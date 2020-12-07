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
// �}�N��
//-----------------------------------------------------------------------------------
//�p���V���[�g�֘A

#define RHO 1.21							//��C���x(1.2kg/m^3)
#define G_ACCE 9.8							//�d�͉����x
#define PI 3.1415							//�~����

#define KITE_RHO 0.2						//���̖��x? ��{�I�ɂ͈����₷���悤200g/1m^2

//�������̎��ʊ֘A
#define EPSILON 0.82
#define K_11 0.5
#define K_33 1.0
#define K_44 0.24
#define K_15 0.75
#define K_66 0

#define SL_NUM 4		//�T�X�y���V�������C���̐�
#define RISER_NUM 4		//���C�U�[�̐�

//�΂ˊ֘A
#define K_SPRING 1000.0		//�΂˒萔
#define D_SPRING 1.0		//�����W��


//---------------------------------------------------------------------------------------------------------------------
//�p���V���[�g�V�~�����[�^
//---------------------------------------------------------------------------------------------------------------------
namespace parachute3d
{

	//! structure "canopy"
	struct canopy
	{
		double m;		//����
		double z;		//������L���m�s�[�̒��S���W
		double R_0;		//�c�����Ă��Ȃ��Ƃ��̔��a
		double R_P;		//�c���������̔��a, R_P = 2/3 R_0
		double S_0;		//���e�ʐ�
		double I_aa;		//�������[�����g
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//�����̏�̒��_(1�̂�)
		Vector3d pos;	//�������\������_(suspension_line�̐������쐬)
		Vector3d local_vel;		//�p���V���[�g���W�n�ł̑��x
		Vector3d local_vel_pre;		//�p���V���[�g���W�n�ł�1step�O�̑��x
		Vector3d global_vel;		//�O���[�o�����W�n�ł̑��x
		Vector3d omega;		//�p���x
		Vector3d omega_pre;		//1step�O�̊p���x
		Vector3d euler_angles;		//�I�C���[�p(�p��)
		Vector3d euler_angles_pre;		//1step�O�̃I�C���[�p
		Vector3d frc;
		Vector3d mom;
	};


	//! structure "suspension_line"
	struct suspension_line
	{
		double m;		//����
		double l;		//�T�X�y���V�������C���̒���
		double z;		//�T�X�y���V�������C���̒��S���W
		double I_aa;		//�������[�����g
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//�L���m�s�[�ƃT�X�y���V�������C���̐ړ_
		Vector3d bottom;	//top���璷��l, riser�Ƃ̐ړ_
		Vector3d pos;
		Vector3d euler_angles;		//�I�C���[�p(�p��)
		Vector3d euler_angles_pre;		//1step�O�̃I�C���[�p
	};

	//! structure "riser"
	struct riser
	{
		double m;		//����
		double z;		//���C�U�[�̒��S���W
		double l;		//���C�U�[�̒���
		double l_min;		//�ʏ펞�̃��C�U�[�̒���
		double l_max;		//�L�΂������̃��C�U�[�̒���
		double I_aa;		//�������[�����g
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//�T�X�y���V�������C���ƃ��C�U�[�̐ړ_
		Vector3d bottom;	//���C�U�[�ƃy�C���[�h�̐ړ_
		Vector3d pos;
	};


	//! structure "payload"
	struct payload
	{
		double m;		//����
		double z;		//�y�C���[�h�̒��S���W
		double a_pl;		//�����̃y�C���[�h�̐��@(m)
		double I_aa;		//�������[�����g
		double I_bb;
		double I_cc;
		double I_xx;
		double I_yy;
		double I_zz;
		Vector3d top;	//���C�U�[�ƃy�C���[�h�̐ړ_
		Vector3d bottom;	//�����̂��\������c��̒��_
		Vector3d pos;
		Vector3d local_vel;		//�p���V���[�g���W�n�ł̑��x
		Vector3d local_vel_pre;		//�p���V���[�g���W�n�ł�1step�O�̑��x
		Vector3d global_vel;		//�O���[�o�����W�n�ł̑��x
		Vector3d global_vel_pre;
		Vector3d omega;		//�p���x
		Vector3d omega_pre;		//1step�O�̊p���x
		Vector3d euler_angles;		//�I�C���[�p(�p��)
		Vector3d euler_angles_pre;		//1step�O�̃I�C���[�p
		Vector3d frc;
		Vector3d mom;
	};

	//! structure "parachute"
	struct parachute
	{
		double m;		//�p���V���[�g�S�̂̎���
		double z_G;		//�p���V���[�g�S�̂̒��S���W
		double z_P;		//���͒��S, z_P = -3/8 * epsilon * R_P
		double I_xx;		//�������[�����g
		double I_yy;
		double I_zz;
		double I_cc;
		Vector3d local_vel;		//�p���V���[�g���W�n�ł̑��x
		Vector3d local_vel_pre;		//�p���V���[�g���W�n�ł�1step�O�̑��x
		Vector3d global_vel;		//�O���[�o�����W�n�ł̑��x
		Vector3d omega;		//�p���x
		Vector3d omega_pre;		//1step�O�̊p���x
		Vector3d euler_angles;		//�I�C���[�p(�p��)
		Vector3d euler_angles_pre;		//1step�O�̃I�C���[�p

		Vector3d pos;		//�O���[�o�����W�n�ł̈ʒu

		Vector3d frc;	//��F
		Vector3d mom;	//���[�����gM
	};


	//������
	void initialize_sim(void);		//�V�~�����[�V�����̏�����
	void initialize_parachute(void);	//�p���V���[�g�p�����[�^�̏�����

	//�p���V���[�g���f���E�ӂ̌v�Z
	void calc_model(double dt);
	//��C�͂̌v�Z
	Vector3d calc_F_aero(void);
	//��̓��[�����g�̌v�Z
	Vector3d calc_M_aero(void);
	//���͂̌v�Z
	Vector3d calc_F_tension(void);
	//���̓��[�����g�̌v�Z
	Vector3d calc_M_tension_can(Vector3d F_tension);
	Vector3d calc_M_tension_pl(Vector3d F_tension);
	//�T�X�y���V�������C���̈ʒu
	void sl_pos(void);

	//��
	void set_wind(double dt);
	//�X�e�b�v��i�߂�
	void step_simulation(double dt);
	//�`��
	void draw_parachute(void);
	//�`��I�v�V����
	void draw_options_1(void);
	void draw_options_2(void);

	void DrawSolidCuboid(Vector3d len);
	void DrawHalfSphere(double rad);
}




#endif

