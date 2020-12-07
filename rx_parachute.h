/*!
 @file rx_parachute.h

 @brief Parachute Simulation

 @author Akari Yabana
 @date 2020
*/
#ifndef _PARACHUTE
#define _PARACHUTE


#include <vector>
#include <Eigen/Dense>
#include <iostream>

//include "rx_utility.h"		// Vector classes
//include "rx_matrix.h"		// Matrix classes
//include "rx_quaternion.h"	// Quaternion classes

#include "rx_solver.h"

using namespace std;
using namespace Eigen;

//-----------------------------------------------------------------------------------
// マクロ
//-----------------------------------------------------------------------------------
//パラシュート関連

#define RHO 1.21							//空気密度(1.2kg/m^3)
#define G_ACCE 9.8							//重力加速度
#define PI 3.1415							//円周率

#define KITE_RHO 0.2						//凧の密度? 基本的には扱いやすいよう200g/1m^2

//見かけの質量関連
#define EPSILON 0.82
#define K_11 0.5
#define K_33 1.0
#define K_44 0.24
#define K_15 0.75
#define K_66 0

#define SL_NUM 4		//サスペンションラインの数
#define RISER_NUM 4		//ライザーの数

//ばね関連
#define K_SPRING 1000.0		//ばね定数
#define D_SPRING 1.0		//減衰係数


//---------------------------------------------------------------------------------------------------------------------
//パラシュートシミュレータ
//---------------------------------------------------------------------------------------------------------------------
namespace parachute3d
{

	//! structure "canopy"
	struct canopy
	{
		double m;		//質量
		double z;		//半球状キャノピーの中心座標
		double R_0;		//膨張していないときの半径
		double R_P;		//膨張した時の半径, R_P = 2/3 R_0
		double S_0;		//投影面積
		double I_aa;		//慣性モーメント
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//屋根の上の頂点(1つのみ)
		Vector3d pos;	//屋根を構成する点(suspension_lineの数だけ作成)
		Vector3d local_vel;		//パラシュート座標系での速度
		Vector3d local_vel_pre;		//パラシュート座標系での1step前の速度
		Vector3d global_vel;		//グローバル座標系での速度
		Vector3d omega;		//角速度
		Vector3d omega_pre;		//1step前の角速度
		Vector3d euler_angles;		//オイラー角(姿勢)
		Vector3d euler_angles_pre;		//1step前のオイラー角
		Vector3d frc;
		Vector3d mom;
	};


	//! structure "suspension_line"
	struct suspension_line
	{
		double m;		//質量
		double l;		//サスペンションラインの長さ
		double z;		//サスペンションラインの中心座標
		double I_aa;		//慣性モーメント
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//キャノピーとサスペンションラインの接点
		Vector3d bottom;	//topから長さl, riserとの接点
		Vector3d pos;
		Vector3d euler_angles;		//オイラー角(姿勢)
		Vector3d euler_angles_pre;		//1step前のオイラー角
	};

	//! structure "riser"
	struct riser
	{
		double m;		//質量
		double z;		//ライザーの中心座標
		double l;		//ライザーの長さ
		double l_min;		//通常時のライザーの長さ
		double l_max;		//伸ばした時のライザーの長さ
		double I_aa;		//慣性モーメント
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//サスペンションラインとライザーの接点
		Vector3d bottom;	//ライザーとペイロードの接点
		Vector3d pos;
	};


	//! structure "payload"
	struct payload
	{
		double m;		//質量
		double z;		//ペイロードの中心座標
		double a_pl;		//立方体ペイロードの寸法(m)
		double I_aa;		//慣性モーメント
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//ライザーとペイロードの接点
		Vector3d bottom;	//立方体を構成する残りの頂点
		Vector3d pos;
		Vector3d local_vel;		//パラシュート座標系での速度
		Vector3d local_vel_pre;		//パラシュート座標系での1step前の速度
		Vector3d global_vel;		//グローバル座標系での速度
		Vector3d global_vel_pre;
		Vector3d omega;		//角速度
		Vector3d omega_pre;		//1step前の角速度
		Vector3d euler_angles;		//オイラー角(姿勢)
		Vector3d euler_angles_pre;		//1step前のオイラー角
		Vector3d frc;
		Vector3d mom;
	};

	//! structure "parachute"
	struct parachute
	{
		double m;		//パラシュート全体の質量
		double z_G;		//パラシュート全体の中心座標
		double z_P;		//圧力中心, z_P = -3/8 * epsilon * R_P
		double I_xx;		//慣性モーメント
		double I_yy;
		double I_zz;
		double I_cc;
		Vector3d local_vel;		//パラシュート座標系での速度
		Vector3d local_vel_pre;		//パラシュート座標系での1step前の速度
		Vector3d global_vel;		//グローバル座標系での速度
		Vector3d omega;		//角速度
		Vector3d omega_pre;		//1step前の角速度
		Vector3d euler_angles;		//オイラー角(姿勢)
		Vector3d euler_angles_pre;		//1step前のオイラー角

		Vector3d pos;		//グローバル座標系での位置

		Vector3d frc;	//力F
		Vector3d mom;	//モーメントM
	};


	//初期化
	void initialize_sim(void);		//シミュレーションの初期化
	void initialize_parachute(void);	//パラシュートパラメータの初期化

	//パラシュートモデル右辺の計算
	void calc_model(double dt);
	//空気力の計算
	Vector3d calc_F_aero(void);
	//空力モーメントの計算
	Vector3d calc_M_aero(void);
	//張力の計算
	Vector3d calc_F_tension(void);
	//張力モーメントの計算
	Vector3d calc_M_tension_can(Vector3d F_tension);
	Vector3d calc_M_tension_pl(Vector3d F_tension);
	//サスペンションラインの位置
	void sl_pos(void);

	//風
	void set_wind(double dt);
	//ステップを進める
	void step_simulation(double dt);
	//描画
	void draw_parachute(void);
	//描画オプション
	void draw_options_1(void);
	void draw_options_2(void);

	void DrawSolidCuboid(Vector3d len);
	void DrawHalfSphere(double rad);
}




#endif

