/*!
 @file parachute.cpp

 @brief Parachute Simulation

 @author Akari Yabana
 @date 2020
*/

// OpenGL
#include <GL/glew.h>
#include <GL/glut.h>

#include "rx_parachute.h"

//軌跡用
double trans_x[5000];
double trans_y[5000];
double trans_z[5000];
double euler_x[5000];
double euler_y[5000];
double euler_z[5000];
int step;

parachute3d::parachute para;		//パラシュート全体の情報格納
parachute3d::canopy can;		//キャノピーの情報格納
parachute3d::suspension_line sl;		//サスペンションラインの情報格納
parachute3d::suspension_line sl_each[SL_NUM];		//サスペンションラインそれぞれの情報格納
parachute3d::riser ris;		//ライザーの情報格納
parachute3d::riser ris_length[RISER_NUM];		//ライザーの長さ情報格納
parachute3d::payload pl;		//ペイロードの情報格納

Vector3d wind;		//風ベクトル
Vector3d wind_a;		//大気流速度ベクトル
Vector3d F_tension;
Vector3d M_tension;

//パラシュートモデルの計算用
//キャノピー右辺
Matrix3d M_m;
Matrix3d Lambda_c;
Vector3d M_m_hat;
Matrix3d I;
Matrix3d H_c;
Vector3d I_hat;
Vector3d M_cr;
Matrix3d rotation_matrix;

Vector3d G;

//ペイロードの右辺
Matrix3d p_M_m;
Matrix3d p_Lambda_c;
Vector3d p_M_m_hat;
Matrix3d p_I;
Matrix3d p_H_c;
Vector3d p_I_hat;
Vector3d p_M_cr;
Matrix3d p_rotation_matrix;

//見かけの質量
double alpha_11;
double alpha_33;
double alpha_44;
double alpha_66;
double alpha_15;

//力
Vector3d F_G;

//モーメント
Vector3d M_G;

//キャノピーとペイロードの距離の初期値
double L;

//オイラー角と角速度の変換行列
Matrix3d euler_omega(double phi, double theta, double psi);
Matrix3d seteuler(double yaw, double pitch, double roll);

using namespace parachute3d;
using namespace std;
using namespace Eigen;


//---------------------------------------------------------------------------------
//パラシュートシミュレータ
//---------------------------------------------------------------------------------

/*!
 * @note パラシュートパラメータの初期化
 */
void parachute3d::initialize_parachute(void)
{
	//canopy-----------------------------
	//can.m = 22.8;		//質量(kg)
	can.m = 1000.8;		//質量(kg)
	can.z = -2.76;		//半球状キャノピーの中心座標(m)
	can.R_0 = 9.75;		//膨張していないときの半径(m)
	can.R_P = 6.50;		//膨張した時の半径(m), R_P = 2/3 R_0
	can.S_0 = PI * pow(can.R_0, 2);		//投影面積
	can.I_aa = 365.93;		//慣性モーメント
	can.I_bb = can.I_aa;
	can.I_cc = 623.98;
	can.I_xx = 539.16;
	can.I_yy = can.I_xx;
	can.I_zz = can.I_cc;
	//------------------------------------

	//suspension line---------------------
	sl.m = 35.3;
	sl.z = 7.50;		//サスペンションラインの中心座標
	sl.I_aa = 1047.09;		//慣性モーメント
	sl.I_bb = sl.I_aa;
	sl.I_cc = 770.76;
	sl.I_xx = 3032.22;
	sl.I_yy = sl.I_xx;
	sl.I_zz = sl.I_cc;
	for (int i = 0; i < SL_NUM; i++) {
		sl_each[i].l = 15.55;		//サスペンションラインの長さ
	}
	//------------------------------------

	//riser-------------------------------
	ris.m = 13.4;
	ris.z = 17.80;		//ライザーの中心座標
	ris.l_min = 5.80;		//通常時のライザーの長さ
	ris.l_max = 7.62;		//伸ばした時のライザーの長さ
	ris.I_aa = 53.43;		//慣性モーメント
	ris.I_bb = ris.I_aa;
	ris.I_cc = 36.98;
	ris.I_xx = 4296.81;
	ris.I_yy = ris.I_xx;
	ris.I_zz = ris.I_cc;
	for (int i = 0; i < RISER_NUM; i++) {
		ris_length[i].l = ris.l_min;
	}
	ris_length[1].l = ris.l_max;
	//-------------------------------------

	//payload------------------------------
	pl.m = 990.0;
	pl.z = 21.20;		//ペイロードの中心座標
	pl.a_pl = 5.22;		//立方体ペイロードの寸法(m) もとは1.22
	pl.I_aa = 245.59;		//慣性モーメント
	pl.I_bb = pl.I_aa;
	pl.I_cc = 245.59;
	pl.I_xx = 445287.93;
	pl.I_yy = pl.I_xx;
	pl.I_zz = pl.I_cc;
	//-------------------------------------

	//parachute----------------------------
	para.m = 1061.5;
	para.z_G = 20.19;		//パラシュート全体の中心座標
	para.z_P = -2.00;		//圧力中心, z_P = -3/8 * epsilon * R_P
	para.I_xx = 453156.13;		//慣性モーメント
	para.I_yy = para.I_xx;
	para.I_zz = 1677.30;
	para.I_cc = para.I_zz;
	//--------------------------------------


	//初期姿勢
	double phi = 0.0;
	double theta =0.5;
	double psi = 0.0;
	para.euler_angles << phi, theta, psi;
	para.euler_angles_pre = para.euler_angles;
	can.euler_angles << phi, theta, psi;
	can.euler_angles_pre = can.euler_angles;
	//ペイロードの初期姿勢
	pl.euler_angles << phi, theta, psi + 0.25 * PI;
	pl.euler_angles_pre = pl.euler_angles;

	//初速
	para.local_vel.setZero();
	para.local_vel_pre.setZero();
	para.global_vel.setZero();
	para.omega.setZero();
	para.omega_pre.setZero();

	can.local_vel.setZero();
	can.local_vel_pre.setZero();
	can.global_vel.setZero();
	can.omega.setZero();
	can.omega_pre.setZero();

	pl.local_vel.setZero();
	pl.local_vel_pre.setZero();
	pl.global_vel.setZero();
	pl.global_vel_pre.setZero();
	pl.omega.setZero();
	pl.omega_pre.setZero();
	
	//位置
	para.pos << 0.0, 0.0, para.z_G;
	can.pos << 0.0, 0.0, can.z;
	sl.pos << 0.0, 0.0, sl.z;
	ris.pos << 0.0, 0.0, ris.z;
	pl.pos << 10.0, 10.0, pl.z-5.0;

	Vector3d tmp0, tmp1, tmp2, tmp3;
	//サスペンションラインの位置　キャノピー側
	tmp0 << 0, can.R_0, 0;
	tmp1 << -can.R_0, 0, 0;
	tmp2 << 0, -can.R_0, 0;
	tmp3 << can.R_0, 0, 0;
	sl_each[0].top = can.pos + tmp0;
	sl_each[1].top = can.pos + tmp1;
	sl_each[2].top = can.pos + tmp2;
	sl_each[3].top = can.pos + tmp3;
	//サスペンションラインの位置　ペイロード側
	tmp0 << 0, 0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
	tmp1 << -0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
	tmp2 << 0, -0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
	tmp3 << 0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
	sl_each[0].bottom = pl.pos + tmp0;
	sl_each[1].bottom = pl.pos + tmp1;
	sl_each[2].bottom = pl.pos + tmp2;
	sl_each[3].bottom = pl.pos + tmp3;
	L = (can.pos - pl.pos).norm();

}
/*!
 * @note その他シミュレーションの初期化
 */
void parachute3d::initialize_sim(void)
{

	//パラシュートモデル右辺の項の初期化
	M_m.setZero();
	Lambda_c.setZero();
	I.setZero();
	H_c.setZero();
	M_m_hat.setZero();
	I_hat.setZero();
	M_cr.setZero();

	p_M_m.setZero();
	p_Lambda_c.setZero();
	p_I.setZero();
	p_H_c.setZero();
	p_M_m_hat.setZero();
	p_I_hat.setZero();
	p_M_cr.setZero();

	//力、モーメントの初期化
	can.frc.setZero();
	can.mom.setZero();
	pl.frc.setZero();
	pl.mom.setZero();

	//見かけの質量alphaの計算
	double m_a = 0.5 * (4 / 3) * PI * RHO * pow(can.R_P, 3) * EPSILON;
	double I_xx_air = (1 / 5) * m_a * pow(can.R_P, 2) * (1 + pow(EPSILON, 2));
	double I_zz_air = (2 / 5) * m_a * pow(can.R_P, 2);

	alpha_11 = K_11 * m_a;
	alpha_33 = K_33 * m_a;
	alpha_44 = K_44 * I_xx_air;
	alpha_66 = K_66 * I_zz_air;
	alpha_15 = K_15 * m_a * para.z_P;

	//グローバル系での重力ベクトル
	G << 0.0, 0.0, G_ACCE;

	//軌跡用配列初期化
	for (int i = 0; i < 1000; i++) {
		trans_x[i] = 0;
		trans_y[i] = 0;
		trans_z[i] = 0;
		euler_x[i] = 0;
		euler_y[i] = 0;
		euler_z[i] = 0;
	}
	step = 0;
}

/*!
 * @note パラシュートモデル右辺の計算
 */
void parachute3d::calc_model(double dt)
{

	//パラシュート全体

	M_m(0,0) = para.m + alpha_11; //対角成分に値を格納
	M_m(1, 1) = para.m + alpha_11;
	M_m(2, 2) = para.m + alpha_33;


	Lambda_c(0,0) = -para.local_vel[1] * para.omega[2];
	Lambda_c(0,1) = para.local_vel[2] * para.omega[1];
	Lambda_c(0,2) = ((para.omega[1] - para.omega_pre[1])/dt) + para.omega[2] * para.omega[0];
	Lambda_c(1,0) = para.local_vel[0] * para.omega[2];
	Lambda_c(1,1) = -para.local_vel[2] * para.omega[0];
	Lambda_c(1,2) = para.omega[1] * para.omega[2] - ((para.omega[0] - para.omega_pre[0])/dt);
	Lambda_c(2,0) = para.local_vel[1] * para.omega[0] - para.local_vel[0] * para.omega[1];
	Lambda_c(2,1) = 0;
	Lambda_c(2,2) = -(pow(para.omega[0], 2) + pow(para.omega[1], 2));
	
	M_m_hat[0] = para.m + alpha_11;
	M_m_hat[1] = para.m + alpha_33;
	M_m_hat[2] = para.m * para.z_G + alpha_15;

	I(0, 0) = para.I_xx + alpha_44;
	I(1, 1) = para.I_yy + alpha_44;
	I(2, 2) = para.I_zz + alpha_66; //対角成分に値を格納

	H_c(0,0) = para.local_vel[2] * para.omega[0] - para.local_vel[0] * para.omega[2] - ((para.local_vel[1] - para.local_vel_pre[1]) / dt);
	H_c(0,1) = -para.omega[1] * para.omega[2];
	H_c(0,2) = para.local_vel[1] * para.local_vel[2];
	H_c(1,0) = para.local_vel[2] * para.omega[1] - para.local_vel[1] * para.omega[2] + ((para.local_vel[0] - para.local_vel_pre[0]) / dt);
	H_c(1,1) = para.omega[0] * para.omega[2];
	H_c(1,2) = -para.local_vel[0] * para.local_vel[2];

	I_hat[0] = para.m * para.z_G + alpha_15;
	I_hat[1] = para.I_yy + alpha_44 - para.I_zz - alpha_66;
	I_hat[2] = alpha_33 - alpha_11;

	M_cr[2] = (para.I_yy - para.I_xx) * para.omega[0] * para.omega[1];

	rotation_matrix = seteuler(para.euler_angles[0], para.euler_angles[1], para.euler_angles[2]);
}

/*!
 * @note 空気力の計算
 */
Vector3d parachute3d::calc_F_aero(void)
{
	Vector3d wind_b = -rotation_matrix.transpose() * wind; //風ベクトルを物体座標系にコンバート
	Vector3d Va = -can.local_vel - wind_b; //対気速度ベクトル

	double alpha_sp = acos(Va[2] / (sqrt(pow(Va[0], 2) + pow(Va[1], 2) + pow(Va[2], 2))));
	alpha_sp = alpha_sp * (180 / PI);
	double CD = calc_CD(alpha_sp);
	double q = 0.5 * RHO * Va.norm();

	Vector3d F_aero = CD * q * can.S_0 * (Va / Va.norm());
	return F_aero;
}


/*!
 * @note 空力モーメントの計算
 */
Vector3d parachute3d::calc_M_aero(void)
{
	Vector3d wind_b = -rotation_matrix.transpose() * wind; //風ベクトルを物体座標系にコンバート
	Vector3d Va = can.local_vel - wind_b; //対気速度ベクトル

	double alpha = atan(Va[0] / Va[2]);
	alpha = alpha * (180 / PI);
	double beta = atan(Va[2] / sqrt(pow(Va[0], 2) + pow(Va[2], 2)));
	beta = beta * (180 / PI);

	double Cn;	//ライザーの長さがすべて等しいときCn=0. 長さが違うとき、パラシュートが回転する

	double lk_0 = (ris_length[0].l - ris.l_min) / (ris.l_max - ris.l_min);
	double lk_1 = (ris_length[1].l - ris.l_min) / (ris.l_max - ris.l_min);
	Cn = 0;
	Vector3d C;
	C << calc_Cm(beta), calc_Cm(alpha), Cn; //ライザーの長さがすべて等しいときCn=0

	double q = 0.5 * RHO * Va.norm();
	Vector3d M_aero = 2 * q * can.S_0 * can.R_0 * C;

	return M_aero;
}

Vector3d parachute3d::calc_F_tension(void)
{
	Vector3d dx = can.pos - pl.pos; //キャノピーとペイロードの距離
	Vector3d vel = can.global_vel - pl.global_vel; //キャノピーとペイロードの速度差

	//張力計算
	F_tension = K_SPRING * (dx.norm() - L) * (dx / dx.norm());
	return F_tension;
}

Vector3d parachute3d::calc_M_tension_can(Vector3d F_tension)
{
	Vector3d M_tension_can, candir[4];

		candir[0] << pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		candir[1] << pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		candir[2] << -pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		candir[3] << -pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		for (int i = 0; i < 4; i++) {
			M_tension_can += -F_tension.cross(candir[i]);
		}

		return M_tension_can;
}

Vector3d parachute3d::calc_M_tension_pl(Vector3d F_tension)
{
	Vector3d M_tension_pl, pldir[4];

		pldir[0] << pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		pldir[1] << pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		pldir[2] << -pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		pldir[3] << -pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;

		for (int i = 0; i < 4; i++) {
			M_tension_pl += F_tension.cross(pldir[i]);
		}
		return M_tension_pl;
}

void parachute3d::sl_pos(void)
{
		Vector3d tmp0, tmp1, tmp2, tmp3;
		tmp0 << 0, can.R_0, 0;
		tmp1 << -can.R_0, 0, 0;
		tmp2 << 0, -can.R_0, 0;
		tmp3 << can.R_0, 0, 0;
		sl_each[0].top = can.mom * (tmp0 + can.pos);
		sl_each[1].top = can.mom * (tmp1 + can.pos);
		sl_each[2].top = can.mom * (tmp2 + can.pos);
		sl_each[3].top = can.mom * (tmp3 + can.pos);



		tmp0 << 0, 0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
		tmp1 << -0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
		tmp2 << 0, -0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
		tmp3 << 0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
		sl_each[0].bottom = pl.mom * (tmp0 + pl.pos);
		sl_each[1].bottom = pl.mom * (tmp1 + pl.pos);
		sl_each[2].bottom = pl.mom * (tmp2 + pl.pos);
		sl_each[3].bottom = pl.mom * (tmp3 + pl.pos);

		return;
}

void
parachute3d::set_wind(double dt)
{
	wind << 0.00005, 0.0, 0.0;
}


/*!
 * @note シミュレーションを進める
 */
void
parachute3d::step_simulation(double dt)
{
	parachute3d::set_wind(dt);//風のセット

	//モデルの計算
	parachute3d::calc_model(dt);


	//--------------- キャノピー --------------------------

	Vector3d F_aero = parachute3d::calc_F_aero(); //空気力
	Vector3d M_aero = parachute3d::calc_M_aero(); //空力モーメント
	Vector3d can_F_G = can.m * G; //重力（グローバル系）
	can_F_G = rotation_matrix.transpose() * can_F_G; //重力をローカル系にコンバート

	can.frc = F_aero + can_F_G;
	can.mom += M_aero;

	//速度
	can.local_vel = can.local_vel_pre + dt * M_m.inverse() * (can.frc - Lambda_c * M_m_hat);
	can.omega = can.omega_pre + dt * I.inverse() * (can.mom - H_c * I_hat - M_cr);
	can.global_vel = rotation_matrix * can.local_vel;
	
	//姿勢
	double phi = can.euler_angles[0];
	double theta = can.euler_angles[1];
	double psi = can.euler_angles[2];
	Matrix3d e_o;
	e_o = euler_omega(phi, theta, psi);
	can.euler_angles = can.euler_angles_pre + dt * (e_o * can.omega); //オイラー角
	//位置
	can.pos += can.global_vel*dt + 0.5*((F_aero / can.m) + G)*dt*dt;


	//速度、姿勢の更新
	//can.global_vel += ((F_aero / can.m) + G)*dt;	
	// can.local_vel_pre = can.local_vel;
	// can.euler_angles_pre = can.euler_angles;

	//--------------- ペイロード -----------------------------

	Vector3d pl_F_G = pl.m * G; //重力（グローバル系）
	//pl_F_G = p_rotation_matrix.transpose() * pl_F_G; //重力をローカル系にコンバート

	pl.frc = pl_F_G;

	//速度
	pl.global_vel = pl.global_vel_pre + dt * (pl_F_G / pl.m);
	//pl.omega += pl.mom.inverse() * (kite.mom - cross(kite.omega, I_omega) - Vec3(DAMP*kite.omega.data[0], DAMP*kite.omega.data[1], DAMP_yaw*kite.omega.data[2]))*dt;
	//姿勢
	double p_phi = pl.euler_angles[0];
	double p_theta = pl.euler_angles[1];
	double p_psi = pl.euler_angles[2];
	Matrix3d p_e_o;
	p_e_o = euler_omega(p_phi, p_theta, p_psi);
	pl.euler_angles = pl.euler_angles_pre + dt * (p_e_o * pl.omega); //オイラー角
	//位置
	pl.pos += pl.global_vel*dt + 0.5*(pl_F_G / pl.m)*dt*dt;


//----------- 引張力 ------------------
		// Vector3d dx = can.pos - pl.pos; //キャノピーとペイロードの距離
		// Vector3d vel = can.global_vel - pl.global_vel; //キャノピーとペイロードの速度差

		// //張力計算
		// F_tension = K_SPRING * (dx.norm() - L) * (dx / dx.norm());
	Vector3d F_tension_can, F_tension_pl;
		F_tension_can = -calc_F_tension();
		F_tension_pl = calc_F_tension();

		//張力によるモーメント
		// Vector3d M_tension_can, M_tension_pl, candir[4], pldir[4];
		// pldir[0] << pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		// pldir[1] << pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		// pldir[2] << -pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		// pldir[3] << -pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;

		// candir[0] << pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		// candir[1] << pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		// candir[2] << -pl.a_pl / 2, pl.a_pl / 2, pl.a_pl / 2;
		// candir[3] << -pl.a_pl / 2, -pl.a_pl / 2, pl.a_pl / 2;
		// for (int i = 0; i < 4; i++) {
		// 	M_tension_pl += F_tension.cross(pldir[i]);
		// 	M_tension_can += -F_tension.cross(candir[i]);
		// }
		Vector3d M_tension_can = calc_M_tension_can(F_tension_can);
		Vector3d M_tension_pl = calc_M_tension_pl(F_tension_pl);


		//cout << "mae:" << dx.norm() - L << endl;
		

		//位置、速度、姿勢の修正
		Vector3d vel_sp_can = dt * F_tension_can;
		Vector3d vel_sp_pl = dt * F_tension_pl;
		can.frc += F_tension_can;
		pl.frc += F_tension_pl;
		can.pos += vel_sp_can * dt;
		pl.pos += vel_sp_pl * dt;
		can.global_vel += vel_sp_can;
		can.local_vel += rotation_matrix.transpose() * vel_sp_can;
		pl.global_vel += vel_sp_pl;
		can.mom += M_tension_can;
		pl.mom += M_tension_pl;

		Vector3d can_I_omega;
		can_I_omega = I * can.omega;
		can.omega+=I.inverse()*(can.mom-can.omega.cross(can_I_omega))*dt;
		Matrix3d pl_I;
		pl_I(0, 0) = (pl.m * pow(pl.a_pl, 2)) / 6;
		pl_I(1, 1) = (pl.m * pow(pl.a_pl, 2)) / 6;
		pl_I(2, 2) = (pl.m * pow(pl.a_pl, 2)) / 6;
		Vector3d pl_I_omega = pl_I * pl.omega;
		pl.omega+=pl_I.inverse()*(pl.mom-pl.omega.cross(pl_I_omega))*dt;

		//姿勢
		Matrix3d e_o;
		e_o = euler_omega(can.euler_angles[0], can.euler_angles[1], can.euler_angles[2]);
		can.euler_angles = can.euler_angles_pre + dt * (e_o * can.omega); //オイラー角
		Matrix3d pl_e_o;
		e_o = euler_omega(pl.euler_angles[0], pl.euler_angles[1], pl.euler_angles[2]);
		pl.euler_angles = pl.euler_angles_pre + dt * (pl_e_o * pl.omega); //オイラー角
		//cout << "ato:" <<  dx2.norm() - L << endl;

		sl_pos();

		// Vector3d tmp0, tmp1, tmp2, tmp3;
		// tmp0 << 0, can.R_0, 0;
		// tmp1 << -can.R_0, 0, 0;
		// tmp2 << 0, -can.R_0, 0;
		// tmp3 << can.R_0, 0, 0;
		// sl_each[0].top = can.pos + tmp0;
		// sl_each[1].top = can.pos + tmp1;
		// sl_each[2].top = can.pos + tmp2;
		// sl_each[3].top = can.pos + tmp3;

		// tmp0 << 0, 0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
		// tmp1 << -0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
		// tmp2 << 0, -0.5 * pl.a_pl * sqrt(2.0), -0.5 * pl.a_pl;
		// tmp3 << 0.5 * pl.a_pl * sqrt(2.0), 0, -0.5 * pl.a_pl;
		// sl_each[0].bottom = pl.pos + tmp0;
		// sl_each[1].bottom = pl.pos + tmp1;
		// sl_each[2].bottom = pl.pos + tmp2;
		// sl_each[3].bottom = pl.pos + tmp3;



	//速度、姿勢の更新
	pl.global_vel_pre = pl.global_vel;
	can.local_vel_pre = can.local_vel;
	pl.euler_angles_pre = pl.euler_angles;



	//cout << para.local_vel.data[0] << "," << para.local_vel.data[1] << "," << para.local_vel.data[2] << endl;
	//para.local_vel = para.local_vel_pre + dt * M_m.inverse() * (F_aero + F_G - Lambda_c * M_m_hat);
	//cout << para.local_vel.data[0] << "," << para.local_vel.data[1] << "," << para.local_vel.data[2] << endl;

	//para.omega = para.omega_pre + dt * I.inverse() * (M_aero - H_c * I_hat - M_cr);

	//para.global_vel = rotation_matrix * para.local_vel;

	//double phi = para.euler_angles[0];
	//double theta = para.euler_angles[1];
	//double psi = para.euler_angles[2];
	////cout << "stepsim: " << phi << "," << theta << "," << psi << endl;
	//Matrix3d e_o;
	//e_o = euler_omega(phi, theta, psi);


	//para.euler_angles = para.euler_angles_pre + dt * (e_o * para.omega); //オイラー角

	//cout << para.pos.data[0] << "," << para.pos.data[1] << "," << para.pos.data[2] << endl;
	//cout << para.global_vel.data[0] << "," << para.global_vel.data[1] << "," << para.global_vel.data[2] << endl;

	//para.pos += para.global_vel*dt + 0.5*((F_canopy / para.m) + g)*dt*dt;	//位置の更新
	//can.pos += para.global_vel*dt + 0.5*para.frc*dt*dt / para.m;	//位置の更新
	//sl.pos += para.global_vel*dt + 0.5*para.frc*dt*dt / para.m;	//位置の更新
	//ris.pos += para.global_vel*dt + 0.5*para.frc*dt*dt / para.m;	//位置の更新
	//pl.pos += para.global_vel*dt + 0.5 * (wind + g) *dt*dt;
	//para.global_vel += ((F_canopy / para.m) + g)*dt;						//速度の更新
	//cout << para.global_vel << endl;
	//para.local_vel_pre = para.local_vel;
	//para.euler_angles_pre = para.euler_angles;
	/*if (ris_length[1].l-0.003 > ris.l_min) {
		ris_length[1].l -= 0.003;
	}
	else {
		ris_length[1].l = ris.l_min;
	}*/

}

/*!
 * @note パラシュート描画
 */
void
parachute3d::draw_parachute(void)
{


	double scale = 0.2;

	//radianからdegreeに変換
	double rotate_angles[3];
	double p_rotate_angles[3];
	double s_rotate_angles[3];
	for (int i = 0; i < 3; i++) {
		rotate_angles[i] = can.euler_angles[i] * 180 / PI;
		p_rotate_angles[i] = pl.euler_angles[i] * 180 / PI;
		s_rotate_angles[i] = sl.euler_angles[i] * 180 / PI;
	}



	//cout << "degree:" << rotate_angles[0] << "," << rotate_angles[1] << "," << rotate_angles[2] << endl;
	//cout << "radian:" << para.euler_angles[0] << "," << para.euler_angles[1] << "," << para.euler_angles[2] << endl;

	//canopy
	glPushMatrix();

	glColor3f(1.0f, 1.0f, 1.0f);//white
	
	glTranslated(can.pos[0] * scale, -can.pos[2] * scale, can.pos[1] * scale);
	//glTranslated(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	glRotated(rotate_angles[0], 1.0, 0.0, 0.0);
	glRotated(-rotate_angles[2], 0.0, 1.0, 0.0);
	glRotated(rotate_angles[1], 0.0, 0.0, 1.0);

	DrawHalfSphere(can.R_0 * scale);

	//glColor3f(0.5f, 0.5f, 0.5f);//gray
	//glBegin(GL_LINES);
	//glVertex3d(rad * 0.5, 0, 0);
	//glVertex3d(0, -(sl_length[0].l + ris_length[0].l) * scale, 0);
	//glVertex3d(-rad * 0.5, 0, 0);
	//glVertex3d(0, -(sl_length[0].l + ris_length[0].l) * scale, 0);
	//glVertex3d(0, 0, rad * 0.5);
	//glVertex3d(0, -(sl_length[0].l + ris_length[0].l) * scale, 0);
	//glVertex3d(0, 0, -rad * 0.5);
	//glVertex3d(0, -(sl_length[0].l + ris_length[0].l) * scale, 0);
	//glEnd();

	glPopMatrix();

	

	glPushMatrix();
	//payload
	glColor3f(1.0f, 1.0f, 1.0f);//white
	//glTranslated(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	glTranslated(pl.pos[0] * scale, -pl.pos[2] * scale, pl.pos[1] * scale);
	glRotated(p_rotate_angles[0], 1.0, 0.0, 0.0);
	glRotated(-p_rotate_angles[2], 0.0, 1.0, 0.0);
	glRotated(p_rotate_angles[1], 0.0, 0.0, 1.0);
	//cout << para.euler_angles[0] << para.euler_angles[1] << para.euler_angles[2] << endl;
	//glTranslated(pl.pos[0] * scale, -pl.pos[2] * scale, pl.pos[1] * scale);
	glutSolidCube(pl.a_pl * scale);

	//glColor3f(1.0f, 1.0f, 1.0f);//white
	//glTranslated(pl.pos[0] * scale, -pl.pos[2] * scale, pl.pos[1] * scale);
	//glRotated(rotate_angles[0], 1.0, 0.0, 0.0);
	//glRotated(-rotate_angles[2], 0.0, 1.0, 0.0);
	//glRotated(rotate_angles[1], 0.0, 0.0, 1.0);
	////cout << para.euler_angles[0] << para.euler_angles[1] << para.euler_angles[2] << endl;
	//glutSolidCube(1.0);

	glPopMatrix();

	//line
	//glPushMatrix();
	//glColor3f(0.0f, 0.0f, 0.0f);//black
	//glBegin(GL_LINES);
	//glVertex3d(pl.pos[0] * scale, -pl.pos[2] * scale, pl.pos[1] * scale);
	//glVertex3d(can.pos[0] * scale, -can.pos[2] * scale, can.pos[1] * scale);
	//glEnd();
	//glPopMatrix();
	glPushMatrix();
	glColor3f(0.0f, 0.0f, 0.0f);//black
	//glRotated(s_rotate_angles[0], 1.0, 0.0, 0.0);
	//glRotated(-s_rotate_angles[2], 0.0, 1.0, 0.0);
	//glRotated(s_rotate_angles[1], 0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	for (int i = 0; i < SL_NUM; i++) {
		glPushMatrix();
		glVertex3d(sl_each[i].bottom[0], -sl_each[i].bottom[2], sl_each[i].bottom[1]);
		glVertex3d(sl_each[i].top[0], -sl_each[i].top[2], sl_each[i].top[1]);
		glPopMatrix();
	}
	glEnd();
	glPopMatrix();

	//trans_x[step] = para.pos[0] * scale;
	//trans_y[step] = -para.pos[2] * scale;
	//trans_z[step] = para.pos[1] * scale;
	//euler_x[step] = rotate_angles[0];
	//euler_y[step] = rotate_angles[1];
	//euler_z[step] = rotate_angles[2];

	//if (step < 5000) {
	//	glPointSize(1.0);
	//	for (int i = 0; i <= step; i++) {
	//		glPushMatrix();
	//		glColor3f(1.0f, 1.0f, 1.0f);//white
	//		glTranslated(trans_x[i], trans_y[i], trans_z[i]);
	//		glRotated(euler_x[i], 1.0, 0.0, 0.0);
	//		glRotated(-euler_z[i], 0.0, 1.0, 0.0);
	//		glRotated(euler_y[i], 0.0, 0.0, 1.0);
	//		glBegin(GL_POINTS);
	//		glVertex3d(0, 0, 0);
	//		glEnd();
	//		glPopMatrix();
	//	}

	//	step++;
	//}
}


/*!
 * @note debug: 力の大きさを矢印で表示
 */
void
parachute3d::draw_options_1(void)
{
	double scale = 0.1;
	glDisable(GL_LIGHTING);

	glLineWidth(2.0f);

	glPushMatrix();
	glTranslated(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	//glRotated(para.euler_angles[0], 1.0, 0.0, 0.0);
	//glRotated(para.euler_angles[1], 0.0, 1.0, 0.0);
	//glRotated(para.euler_angles[2], 0.0, 0.0, 1.0);

	////重力
	//glColor3f(1.0f, 0.0f, 0.0f);//red

	//glBegin(GL_LINES);
	//glVertex3d(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	//glVertex3d(para.pos[0] * scale + scale * G[0], -para.pos[2] * scale - scale * G[2], para.pos[1] * scale + scale * G[1]);
	//glEnd();

	////F_canopy
	//glColor3f(1.0f, 1.0f, 0.0f);//yellow

	//glBegin(GL_LINES);
	//glVertex3d(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	//glVertex3d(para.pos[0] * scale + scale * F_canopy[0], -para.pos[2] * scale - scale * F_canopy[2], para.pos[1] * scale + scale * F_canopy[1]);
	//glEnd();

	////wind
	//glColor3f(0.0f, 1.0f, 1.0f);
	//glBegin(GL_LINES);
	//glVertex3d(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
	//glVertex3d(para.pos[0] * scale + scale * wind[0], -para.pos[2] * scale - scale * wind[2], para.pos[1] * scale + scale * wind[1]);
	//glEnd();

	//重力
	glColor3f(1.0f, 0.0f, 0.0f);//red

	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(scale * G[0], -scale * G[2], scale * G[1]);
	glEnd();

	//F_canopy
	//glColor3f(1.0f, 1.0f, 0.0f);//yellow

	//glBegin(GL_LINES);
	//glVertex3d(0, 0, 0);
	//glVertex3d(scale * F_aero[0], -scale * F_aero[2], scale * F_aero[1]);
	//glEnd();

	//wind
	glColor3f(0.0f, 1.0f, 1.0f);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(scale * wind[0], -scale * wind[2], scale * wind[1]);
	glEnd();
	glEnable(GL_LIGHTING);

	glPopMatrix();
}
//
//void parachute3d::draw_options_2(void)
//{
//	glPushMatrix();
//	double scale = 0.1;
//	glColor3f(1.0f, 0.0f, 0.0f);
//	glPointSize(10.0);
//	glTranslated(para.pos[0] * scale, -para.pos[2] * scale, para.pos[1] * scale);
//	glRotated(para.euler_angles[0], 1.0, 0.0, 0.0);
//	glRotated(para.euler_angles[1], 0.0, 1.0, 0.0);
//	glRotated(para.euler_angles[2], 0.0, 0.0, 1.0);
//	glTranslatef(pl.pos[0] * scale, -pl.pos[2] * scale, pl.pos[1] * scale);
//
//	glBegin(GL_POINTS);
//	glVertex3d(0, 0, 0);
//	glEnd();
//	glPopMatrix();
//}


/*!
 * @note 角速度からオイラー角への変換行列
 */
Matrix3d euler_omega(double phi, double theta, double psi)
{
	Matrix3d e_o;

	e_o(0,0) = 1;
	e_o(0,1) = sin(phi) * tan(theta);
	e_o(0,2) = cos(phi) * tan(theta);
	e_o(1,0) = 0;
	e_o(1,1) = cos(phi);
	e_o(1,2) = sin(phi);
	e_o(2,0) = 0;
	e_o(2,1) = sin(phi) * (1 / cos(theta));
	e_o(2,2) = cos(phi) * (1/cos(theta));

	return e_o;
}

/*!
 * @note オイラー角から回転行列を設定
 */
Matrix3d seteuler(double yaw, double pitch, double roll)
{
	// yaw is CW around y-axis, pitch is CCW around x-axis, and roll is CW around z-axis
	double cy = cos(yaw);
	double sy = sin(yaw);
	double cp = cos(pitch);
	double sp = sin(pitch);
	double cr = cos(roll);
	double sr = sin(roll);

	double cc = cy * cr;
	double cs = cy * sr;
	double sc = sy * cr;
	double ss = sy * sr;
	Matrix3d euler;
	euler << cc + sp * ss, cs - sp * sc, -sy * cp, -cp * sr, cp*cr, -sp, sc - sp * cs, ss + sp * cc, cy*cp;

	return euler;
}

/*!
 * cenを中心で辺の長さがlenの直方体の描画
 * @param[in] cen 直方体の中心
 * 
 */
void parachute3d::DrawSolidCuboid(Vector3d len)
{
	len *= 0.5;
	Vector3d corner[8];

	corner[0] = Vector3d(-len[0], -len[1], -len[2]);
	corner[1] = Vector3d(-len[0], len[1], -len[2]);
	corner[2] = Vector3d(-len[0], len[1], len[2]);
	corner[3] = Vector3d(-len[0], -len[1], len[2]);

	corner[4] = Vector3d(len[0], -len[1], -len[2]);
	corner[5] = Vector3d(len[0], len[1], -len[2]);
	corner[6] = Vector3d(len[0], len[1], len[2]);
	corner[7] = Vector3d(len[0], -len[1], len[2]);

	int index[6][4] = { { 3, 2, 1, 0 },
						{ 4, 5, 6, 7 },
						{ 3, 0, 4, 7 },
						{ 1, 2, 6, 5 },
						{ 0, 1, 5, 4 },
						{ 2, 3, 7, 6 } };

	glPushMatrix();
	glBegin(GL_QUADS);

	// x軸に垂直な面
	glNormal3d(-1.0, 0.0, 0.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[0][i]][0], corner[index[0][i]][1], corner[index[0][i]][2]);

	glNormal3d(1.0, 0.0, 0.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[1][i]][0], corner[index[1][i]][1], corner[index[1][i]][2]);

	// y軸に垂直な面
	glNormal3d(0.0, -1.0, 0.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[2][i]][0], corner[index[2][i]][1], corner[index[2][i]][2]);
	glNormal3d(0.0, 1.0, 0.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[3][i]][0], corner[index[3][i]][1], corner[index[3][i]][2]);

	// z軸に垂直な面
	glNormal3d(0.0, 0.0, -1.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[4][i]][0], corner[index[4][i]][1], corner[index[4][i]][2]);
	glNormal3d(0.0, 0.0, 1.0);
	for (int i = 0; i < 4; ++i) glVertex3d(corner[index[5][i]][0], corner[index[5][i]][1], corner[index[5][i]][2]);

	glEnd();
	glPopMatrix();
}


/*!
 * 半径radの半球
 * @param[in] rad 半球の半径
 */
void parachute3d::DrawHalfSphere(double rad)
{

	double c[4] = { 0, 1, 0, 0 };
	double c_init[4] = { 0, 0, 0, 0 };
	glEnable(GL_CLIP_PLANE0);

	glPushMatrix();
	glClipPlane(GL_CLIP_PLANE0, c);
	glutSolidSphere(rad, 20, 10);
	glPopMatrix();

	glClipPlane(GL_CLIP_PLANE0, c_init);
	glDisable(GL_CLIP_PLANE0);
}
