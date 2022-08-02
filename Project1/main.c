#include "DllHeader.h"
#include <math.h>

// Inverse Factorial
#define f2          ((double)(1./2.))
#define f3          ((double)(f2/3.))
#define f4          ((double)(f3/4.))
#define f5          ((double)(f4/5.))
#define f6          ((double)(f5/6.))
#define f7          ((double)(f6/7.))
#define f8          ((double)(f7/8.))
#define f9          ((double)(f8/9.))
#define f10         ((double)(f9/10.))
#define f11         ((double)(f10/11.))
#define f12         ((double)(f11/12.))
#define f13         ((double)(f12/13.))
#define f14         ((double)(f13/14.))
#define f15         ((double)(f14/15.))

// Floating-polong Constants
#define PI          ((double)3.1415926535897932384626433832795)
#define SQRT2       ((double)1.414213562)
#define SQRT3       ((double)1.732050808)
#define INV_SQRT3   ((double)0.577350269)
#define INV_SQRT2   ((double)1./SQRT2)
#define INV3        ((double)0.33333333333333333)
#define INV_2PI     ((double)0.15915494309189533576888376337251)

// Macro Functions
#define BOUND_PI(x)	((x>0)?((x)-2.*PI*(long)((x+PI)*INV_2PI)):((x)-2.*PI*(long)((x-PI)*INV_2PI)))
#define SIN(x,x2)	((x)*((double)1.-(x2)*(f3-(x2)*(f5-(x2)*(f7-(x2)*(f9-(x2)*(f11-(x2)*(f13-(x2)*f15))))))))
#define COS(x2)		((double)1.-(x2)*(f2-(x2)*(f4-(x2)*(f6-(x2)*(f8-(x2)*(f10-(x2)*(f12-(x2)*f14)))))))
#define Rm2Rpm      ((double)9.5492965855137201461330258023509)
#define Rpm2Rm      ((double)0.10471975511965977461542144610932)


// Global vairable
double Tsamp = 0.;
double InvTsamp = 0.;
double Vdc = 0.;

//SC
int scCall_cnt = 0, scCall = 10;

// Speed controller variables for PWM0
double Wrpm_ref = 0., Wrpm_ref_set = 0., Wrm_ref = 0., Wr_ref = 0., Wr_enc = 0.;
double Err_Wrm = 0., Te_ref_integ = 0., Kp_sc = 0., Ki_sc = 0., Ki_scT = 0., Ka_sc = 0., Te_ref = 0., Te_real = 0.;
double Te_ref_max = 65.1, Iqse_ref_max = 55.15;
double Kt = 0., Te_Anti = 0., Te_ref_set = 0.;
double Wrm2Wrpm = 30. / PI, Wc_sc = 2.*PI*2.;

// Pole voltage reference
double Van_ref = 0., Vbn_ref = 0., Vcn_ref = 0.;
double Van1_ref = 0., Vbn1_ref = 0., Vcn1_ref = 0.;
double Van_ref_set = 0., Vbn_ref_set = 0., Vcn_ref_set = 0.;

// Anti-wind up voltage
double Vdse_Anti = 0., Vqse_Anti = 0., Vdse_Sat = 0., Vqse_Sat = 0., Vdss_Sat = 0., Vqss_Sat = 0.;
double Vdse1_Anti = 0., Vqse1_Anti = 0., Vdse1_Sat = 0., Vqse1_Sat = 0., Vdss1_Sat = 0., Vqss1_Sat = 0.;

// Current controller variables for PWM0
double Ias = 0., Ibs = 0., Ics = 0.;
double Idss = 0., Iqss = 0., Idse = 0., Iqse = 0., Idse_ref_set = 0., Iqse_ref_set = 0.;
double Idss_ref_set = 0., Iqss_ref_set = 0., Idss_ref = 0., Iqss_ref = 0.;
double Idse_ref = 0., Iqse_ref = 0.;
double Err_Idse = 0., Vdse_ref = 0., Vdss_ref = 0., Vdse_ref_fb = 0., Vdse_ref_ff = 0., Vdse_ref_integ = 0.;
double Err_Iqse = 0., Vqse_ref = 0., Vqss_ref = 0., Vqse_ref_fb = 0., Vqse_ref_ff = 0., Vqse_ref_integ = 0.;
double Vas_ref = 0., Vbs_ref = 0., Vcs_ref = 0.;
double Wc_cc = 2.*PI*250., Kpd_cc = 0., Kid_cc = 0., Kpq_cc = 0., Kiq_cc = 0., Kad_cc = 0., Kaq_cc = 0.;

double Vmax = 0., Vmin = 0., Vsn = 0., Vsn_min = 0., Vsn_max = 0.;
double Cos_Thetar = 0., Sin_Thetar = 0., Theta_set = 0., Thetar = 0., Thetar_err = 0.;

int Flag_adv = 0;
double Wrc = 0.;

int Flag_DTC = 0, Flag_ZCC = 0;

// Parameters
double PolePair = 2., INV_PolePair = 1. / 2.;
//double Rs = 0.014, Ld = 0.000098, Lq = 0.000108, Lavg = 0., Ldelta = 0.;
//double Rs = 0.014, Ld = 0.00294, Lq = 0.00324, Lavg = 0., Ldelta = 0.;
//double LAMpm = 0.00955;
//double Jeq = 0.0001;----

double Rs = 0.1398, Ld = 3.59e-3, Lq = 4.32e-3, LAMpm = 0.2625, Lavg = 0., Ldelta = 0.;
double Jeq = 0.1;

double Thetar_FB = 0., Thetar_enc = 0.;
double Wr_FB = 0., Wrm_FB = 0., Wrpm_FB = 0.;

// 전류제어, 속도제어
int Flag_INV_sc = 2;		// 토크 제어: 1, 속도 제어: 2
int Flag_PM_CC = 0;			// 전류 제어: 1
int Flag_TS_offset = 0;		// 토크센서 오프셋 플래그

double on_off_state = 0;

double Van_comp = 0., Vbn_comp = 0., Vcn_comp = 0., Vdt = 0.;
double Van_ref_out = 0., Vbn_ref_out = 0., Vcn_ref_out = 0.;


// LPF 1
double wcutLPF1 = 0., w_LPF1 = 0.;
double xLF10 = 0., xLF11 = 0., yLF10 = 0., yLF11 = 0.;
double D0LF1 = 0., D1LF1 = 0., N0LF1 = 0., N1LF1 = 0.;


/*DLL input*/
#define in_Ias					aState->inputs[0]
#define in_Ibs					aState->inputs[1]
#define in_Ics					aState->inputs[2]
#define in_Wrpm_ref				aState->inputs[3]
#define in_Wrm_enc				aState->inputs[4]
#define in_Thetarm_enc			aState->inputs[5]
#define in_Vdc					aState->inputs[6]
#define in_CTRDIR				aState->inputs[7]

/*DLL output*/
#define	out_Van_ref				aState->outputs[0]
#define	out_Vbn_ref				aState->outputs[1]
#define	out_Vcn_ref				aState->outputs[2]

#define out_Vdse_Sat			aState->outputs[3]
#define out_Vqse_Sat			aState->outputs[4]

#define out_Vdse_Anti			aState->outputs[5]
#define out_Vqse_Anti			aState->outputs[6]

#define	out_Idss_DLL			aState->outputs[7]
#define	out_Iqss_DLL			aState->outputs[8]
#define	out_Idse_DLL			aState->outputs[9]
#define	out_Iqse_DLL			aState->outputs[10]

#define	out_Idse_ref			aState->outputs[11]
#define	out_Iqse_ref			aState->outputs[12]

#define	out_Te_ref				aState->outputs[13]
#define	out_Te_sat				aState->outputs[14]

#define	out_Thetar_DLL			aState->outputs[15]
#define out_Thetar_FB			aState->outputs[16]
#define	out_Wr_DLL				aState->outputs[17]
#define	out_Wr_FB				aState->outputs[18]
#define	out_Wr_ref_DLL			aState->outputs[19]

#define out_Ias					aState->outputs[20]
#define out_Ibs					aState->outputs[21]
#define out_Ics					aState->outputs[22]

#define out_Vdse_ref_integ		aState->outputs[23]
#define out_Vqse_ref_integ		aState->outputs[24]

DLLEXPORT void plecsSetSizes(struct SimulationSizes* aSizes)
{
	aSizes->numInputs = 8;
	aSizes->numOutputs = 25;
	aSizes->numStates = 1;
	aSizes->numParameters = 1; //number of user parameters passed in
}

//This function is automatically called at the beginning of the simulation
DLLEXPORT void plecsStart(struct SimulationState* aState)
{
	Tsamp = aState->parameters[0];
	InvTsamp = 1. / Tsamp;

	// speed control parameter
	Wc_sc = 2. * PI * 4.;
	Kp_sc = Jeq * Wc_sc;
	Ki_sc = Jeq * (Wc_sc * Wc_sc) / 5.;
	Ka_sc = 1. / Kp_sc;

	//current control parameter
	Wc_cc = 2. * PI * 250.;
	Kpd_cc = Ld * Wc_cc;
	Kid_cc = Rs * Wc_cc;
	Kpq_cc = Lq * Wc_cc;
	Kiq_cc = Rs * Wc_cc;
	Kad_cc = 1. / Kpd_cc;
	Kaq_cc = 1. / Kpq_cc;

	//전류 주입
	Idse_ref_set = 0.;
	Iqse_ref_set = 0.;
/*
	//LPF 예시
	wcutLPF1 = 2. * PI * 20;
	w_LPF1 = 2. / Tsamp * tan(0.5 * wcutLPF1 * Tsamp);
	N0LF1 = Tsamp * w_LPF1 / (2. + Tsamp * w_LPF1);
	N1LF1 = N0LF1;
	D0LF1 = 1.;
	D1LF1 = -(2. - Tsamp * w_LPF1) / (2. + Tsamp * w_LPF1);
*/
	//Example error message box
	//if (aState->parameters[0] < 0 )
	//	aState->errorMessage = "kp cannot be less than 0";
}


//This function is automatically called every sample time
//output is written to DLL output port after the output delay
DLLEXPORT void plecsOutput(struct SimulationState* aState)
{
	// 출력으로 볼 변수들 설정 (DA)
	out_Van_ref = Van_ref;
	out_Vbn_ref = Vbn_ref;
	out_Vcn_ref = Vcn_ref;

	out_Vdse_Sat = Vdse_Sat;
	out_Vqse_Sat = Vqse_Sat;

	out_Vdse_Anti = Vdse_Anti;
	out_Vqse_Anti = Vqse_Anti;

	out_Idss_DLL = Idss;
	out_Iqss_DLL = Iqss;
	out_Idse_DLL = Idse;
	out_Iqse_DLL = Iqse;

	out_Idse_ref = Idse_ref;
	out_Iqse_ref = Iqse_ref;

	out_Te_ref = Te_ref;
	out_Te_sat = Te_real;

	out_Thetar_DLL = Thetar_enc;

	out_Thetar_FB = Thetar_FB;


	out_Wr_DLL = Wr_enc;

	out_Wr_FB = Wr_FB;
	out_Wr_ref_DLL = Wr_ref;

	out_Ias = Ias;
	out_Ibs = Ibs;
	out_Ics = Ics;

	out_Vdse_ref_integ = Vdse_ref_integ;
	out_Vqse_ref_integ = Vqse_ref_integ;

	on_off_state = in_CTRDIR;

	// 입력받을 변수들 설정 (readAngle)
	Thetar_enc = BOUND_PI(in_Thetarm_enc);
	Thetar_FB = BOUND_PI(Thetar_enc*PolePair);

	Wr_enc = in_Wrm_enc * PolePair;
	Wr_FB = in_Wrm_enc * PolePair;
	Wrm_FB = in_Wrm_enc;
	Wrpm_FB = Wrm_FB * Rm2Rpm;

	Wrm_ref = in_Wrpm_ref * Rpm2Rm;
	Wr_ref = Wrm_ref * PolePair;

	// 고정자 좌표계 상의 전류 입력 받아야 함 (ADC)

	Ias = in_Ias;
	Ibs = in_Ibs;
	Ics = in_Ics;

	/* LPF 예시
	xLF10 = Vout_TS;
	yLF10 = (-D1LF1 * yLF11 + N0LF1 * xLF10 + N1LF1 * xLF11) / D0LF1;
	yLF11 = yLF10;
	xLF11 = xLF10;
	*/

	// 고정자 좌표계 상의 전류 입력 받아야 함 (ADC)
	Idss = (2 * in_Ias - in_Ibs - in_Ics) * INV3;
	Iqss = (in_Ibs - in_Ics) * INV_SQRT3;

	Vdc = in_Vdc;

	double a, b, c;

	a = Thetar_FB * Thetar_FB;

	Cos_Thetar = COS(a);
	Sin_Thetar = SIN(Thetar_FB, a);

	Idse = Cos_Thetar * Idss + Sin_Thetar * Iqss;
	Iqse = -Sin_Thetar * Idss + Cos_Thetar * Iqss;

	/////////////////////////////////
	// Current Controller for PWM0 //
	/////////////////////////////////

	//전류제어기
	// D-axis Synchronous PI Controller
	Err_Idse = Idse_ref - Idse;
	Vdse_ref_integ += Kid_cc * (Err_Idse - Kad_cc * Vdse_Anti)*Tsamp;
	Vdse_ref_fb = Vdse_ref_integ + Kpd_cc * Err_Idse;

	// Q-axis Synchronous PI Controller
	Err_Iqse = Iqse_ref - Iqse;
	Vqse_ref_integ += Kiq_cc * (Err_Iqse - Kaq_cc * Vqse_Anti)*Tsamp;
	Vqse_ref_fb = Vqse_ref_integ + Kpq_cc * Err_Iqse;



	// Feed-forward voltages
	 Vdse_ref_ff = Idse*Rs - Wr_FB *Lq*Iqse;
	Vqse_ref_ff = Iqse*Rs + LAMpm* Wr_FB + Wr_FB *Ld*Idse;
	//Vdse_ref_ff = 0.;
	//Vqse_ref_ff = 0.;

	// d-q Voltage References
	Vdse_ref = Vdse_ref_fb + Vdse_ref_ff;
	Vqse_ref = Vqse_ref_fb + Vqse_ref_ff;

	// 1.5Tsamp*Wr voltage angle compensation
	if (Flag_adv)
	{
		Wrc = Wr_FB;
		b = (Thetar_FB + 1.5*Tsamp*Wrc)*(Thetar_FB + 1.5*Tsamp*Wrc);

		Cos_Thetar = COS(b);
		Sin_Thetar = SIN((Thetar_FB + 1.5*Tsamp*Wrc), b);
	}

	Vdss_ref = Cos_Thetar * Vdse_ref - Sin_Thetar * Vqse_ref;
	Vqss_ref = Sin_Thetar * Vdse_ref + Cos_Thetar * Vqse_ref;

	// stationary to three-phase conversion
	Vas_ref = Vdss_ref;
	Vbs_ref = -0.5*(Vdss_ref - SQRT3 * Vqss_ref);
	Vcs_ref = -(Vas_ref + Vbs_ref);


	/////////////////////////
	//  space vector PWM   //
	/////////////////////////

	if (Vas_ref > Vbs_ref)
	{
		Vmax = Vas_ref;
		Vmin = Vbs_ref;
	}
	else
	{
		Vmax = Vbs_ref;
		Vmin = Vas_ref;
	}

	if (Vcs_ref > Vmax)  Vmax = Vcs_ref;
	if (Vcs_ref < Vmin)  Vmin = Vcs_ref;

	Vsn_max = 0.5*Vdc - Vmax;
	Vsn_min = -0.5*Vdc - Vmin;

	if (Vsn_max <= Vsn_min)
	{
		a = 0.5*(Vsn_max + Vsn_min);
		Vsn_max = Vsn_min = a;
	}

	Vsn = 0.5*(Vsn_max + Vsn_min);
	//Vsn = Vsn_min;

	// SVPWM
	Van_ref = Vas_ref + Vsn;
	Vbn_ref = Vbs_ref + Vsn;
	Vcn_ref = Vcs_ref + Vsn;

	Van_ref_out = Van_ref;
	Vbn_ref_out = Vbn_ref;
	Vcn_ref_out = Vcn_ref;


	// Minimum error OVM, voltage limiter
	b = 0.5 * Vdc;
	Van_ref = (Van_ref > b) ? b : ((Van_ref < -b) ? -b : Van_ref);
	Vbn_ref = (Vbn_ref > b) ? b : ((Vbn_ref < -b) ? -b : Vbn_ref);
	Vcn_ref = (Vcn_ref > b) ? b : ((Vcn_ref < -b) ? -b : Vcn_ref);

	// Recalculation of Voltage
	Vdss_Sat = (2.*Van_ref - Vbn_ref - Vcn_ref) / 3;

	Vqss_Sat = INV_SQRT3 * (Vbn_ref - Vcn_ref);

	// Limited d-q voltage
	Vdse_Sat = Cos_Thetar * Vdss_Sat + Sin_Thetar * Vqss_Sat;
	Vqse_Sat = -Sin_Thetar * Vdss_Sat + Cos_Thetar * Vqss_Sat;

	// Anti-wind up
	Vdse_Anti = Vdse_ref - Vdse_Sat;
	Vqse_Anti = Vqse_ref - Vqse_Sat;

	scCall_cnt++;

	// Speed  controller
	if (scCall_cnt >= scCall)
	{
		if (Flag_INV_sc == 1)	// 동기기 토크 제어
		{
			if (Flag_PM_CC == 1)	// 전류 제어
			{
				Idse_ref = Idse_ref_set;
				Iqse_ref = Iqse_ref_set;
			}
			else
			{
				Idse_ref = 0.;
				Iqse_ref = 0.;
			}
		}
		else if (Flag_INV_sc == 2)	// 동기기 속도 제어
		{
			Err_Wrm = Wrm_ref - (Wrm_FB);

			Te_ref_integ += Ki_sc * (Err_Wrm - Ka_sc * (Te_ref - Te_real))*Tsamp * 10;
			Te_ref = Te_ref_integ + Kp_sc * Err_Wrm;

			a = ((Te_ref > Te_ref_max) ? Te_ref_max : (Te_ref < -Te_ref_max) ? -Te_ref_max : Te_ref);

			Kt = 1.5*PolePair*LAMpm;    //for SPM

			Idse_ref = 0.;
			Iqse_ref = a / Kt;

			a = Iqse_ref;
			Iqse_ref = ((a > Iqse_ref_max) ? Iqse_ref_max : (a < -Iqse_ref_max) ? -Iqse_ref_max : a);

			Te_real = Kt * Iqse_ref;
		}

		scCall_cnt = 0;
	}
}
