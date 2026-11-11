// gemc headers
#include "uRwell_strip.h"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

vector<uRwell_strip_found> uRwell_strip::FindStrip(G4ThreeVector xyz, double Edep, uRwellConstants uRwellc, double time, bool isProto, bool isXY, string readout)
{
	vector<uRwell_strip_found> strip_found;
	uRwell_strip_found ThisStrip;

	// Convert deposited energy to the number of primary electrons.
	int N_el = 1e6*Edep/uRwellc.w_i;

	if (N_el == 0) {
		ThisStrip.numberID = -15000;
		ThisStrip.weight   = 1;
		ThisStrip.time     = -1;
		strip_found.push_back(ThisStrip);
		return strip_found;
	}

	N_el = G4Poisson(N_el*uRwellc.gain);

	double strip_coordinate = 0.0;
	double half_size = 0.0;
	double z_real = xyz.z();

	if (isXY) {
		// Small Hall C rectangular prototype (G4Box):
		// direct X-Y readout, with NO stereo-angle rotation.
		if (readout == "x") {
			strip_coordinate = xyz.x();
			half_size = uRwellc.Xhalf_base;
		}
		else if (readout == "y") {
			strip_coordinate = xyz.y();
			half_size = uRwellc.Yhalf;
		}
		else {
			ThisStrip.numberID = -15000;
			ThisStrip.weight   = 1;
			ThisStrip.time     = -1;
			strip_found.push_back(ThisStrip);
			return strip_found;
		}
	}
	else {
		// Big trapezoidal uRwell (G4Trap):
		// preserve the original stereo-strip rotation.
		double x_real = xyz.x()*cos(M_PI*uRwellc.get_stereo_angle()/180)
		              + xyz.y()*sin(M_PI*uRwellc.get_stereo_angle()/180);
		double y_real = xyz.y()*cos(M_PI*uRwellc.get_stereo_angle()/180)
		              - xyz.x()*sin(M_PI*uRwellc.get_stereo_angle()/180);
		(void)x_real;
		strip_coordinate = y_real;
		half_size = uRwellc.Yhalf;
	}

	// Local strip numbering:
	// component = 1 ... number_of_strip for each readout layer.
	double pitch = 2.0*half_size/uRwellc.number_of_strip;
	int stripID = (int) floor((strip_coordinate + half_size)/pitch) + 1;

	// Reject hits outside the active strip range.
	if (stripID < 1 || stripID > uRwellc.number_of_strip) {
		ThisStrip.numberID = -15000;
		ThisStrip.weight   = 1;
		ThisStrip.time     = -1;
		strip_found.push_back(ThisStrip);
		return strip_found;
	}

	// Drift-time contribution.
	double time_dz = fabs(-uRwellc.Zhalf/cm + z_real/cm)/uRwellc.v_drift;

	ThisStrip.numberID = stripID;
	ThisStrip.weight   = N_el;
	ThisStrip.time     = time + time_dz + G4RandGauss::shoot(0., uRwellc.sigma_time);

	strip_found.push_back(ThisStrip);
	return strip_found;
}


double uRwell_strip::Weight_td(int strip, double x, double y, double z, uRwellConstants uRwellc, bool isProto){
	double wght;
	if(Build_strip(strip, uRwellc)){
	 wght=( erf((strip_y+uRwellc.get_strip_width(strip, isProto)/2.-y)/uRwellc.sigma_td/sqrt(2))-erf((strip_y-uRwellc.get_strip_width(strip, isProto)/2.-y)/uRwellc.sigma_td/sqrt(2)))*
			 (erf((strip_x+strip_length/2.-x)/uRwellc.sigma_td/sqrt(2))-erf((strip_x-strip_length/2.-x)/uRwellc.sigma_td/sqrt(2)))/2./2.;
	 if (wght<0) wght=-wght;
	}else{
		wght =-1;
	}
	return wght;
}

bool uRwell_strip::Build_strip(int strip, uRwellConstants uRwellc ){
	
	//strip straight line -> y = mx +c;
	double m = tan(M_PI*uRwellc.get_stereo_angle()/180);
	double c = strip*uRwellc.get_strip_pitch()/cos(M_PI*uRwellc.get_stereo_angle()/180);

   // Trapezoid coordinates
	G4ThreeVector A = {-uRwellc.Xhalf_base, -uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector B =  {uRwellc.Xhalf_base, -uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector C = {-uRwellc.Xhalf_Largebase, uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector D =  {uRwellc.Xhalf_Largebase, uRwellc.Yhalf, uRwellc.Zhalf};
	


	// C-------------D //
	//  -------------  //
	//   -----------   //
	//    A-------B   //
	// Intersection points between strip straight line and Trapezoid straight lines
	
	G4ThreeVector AB_strip = intersectionPoint(m,c,A,B);
	G4ThreeVector BD_strip = intersectionPoint(m,c,B,D);
	G4ThreeVector CD_strip = intersectionPoint(m,c,C,D);
	G4ThreeVector AC_strip = intersectionPoint(m,c,A,C);
	
	vector< G4ThreeVector> strip_points ; // intersection point between strip and the trapezoid sides;
	
	// geometrical characteristic
	double length_strip=0;
	G4ThreeVector first_point;
	G4ThreeVector second_point;
	
	// check if the intersection point is on the segment defined by two points (i.e A and B)
	
	 if(uRwellc.get_strip_kind()=="strip_v"){

		if(pointOnsegment(AC_strip, A, C)) {
			first_point=AC_strip;
			if(pointOnsegment(BD_strip, B, D)) second_point = BD_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
			
		}else if(pointOnsegment(AB_strip, A, B)){
			first_point=AB_strip;
			
			if(pointOnsegment(BD_strip, B, D)) second_point = BD_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
			
		}else{
			return false;
		}
	}
	
	
	
	 if(uRwellc.get_strip_kind()=="strip_u"){

		if(pointOnsegment(BD_strip, B, D)){
			first_point=BD_strip;
			
			if(pointOnsegment(AC_strip, A, C)) second_point = AC_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
		}else if (pointOnsegment(AB_strip, A, B)){
			first_point=AB_strip;
			if(pointOnsegment(AC_strip, A, C)) second_point = AC_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
		}else{
			return false;
		}
		
	}
	
	length_strip = cal_length(first_point, second_point);
	
	strip_length = length_strip;
	strip_endpoint1 = first_point;
	strip_endpoint2 = second_point;


	G4ThreeVector strip_endpoint1_stripFrame = change_of_coordinates(strip_endpoint1, uRwellc );
	G4ThreeVector strip_endpoint2_stripFrame = change_of_coordinates(strip_endpoint2, uRwellc );
	
	strip_y = strip_endpoint1_stripFrame.y();
	strip_x = (strip_endpoint1_stripFrame.x() + strip_endpoint2_stripFrame.x())/2;
	
/*
	cout << "strip: "<<strip<< endl;
	cout << strip_endpoint1.x()<< " "<< strip_endpoint1.y()<<" "<<endl;
	cout << strip_endpoint2.x()<< " "<< strip_endpoint2.y()<<" "<<endl;
	cout <<"done"<<endl;
*/

	return true;
}

G4ThreeVector uRwell_strip::change_of_coordinates( G4ThreeVector A, uRwellConstants uRwellc){
	
	G4ThreeVector XYZ;
	XYZ.setX(A.x()*cos(M_PI*uRwellc.get_stereo_angle()/180) + A.y()*sin(M_PI*uRwellc.get_stereo_angle()/180));
    XYZ.setY(A.y()*cos(M_PI*uRwellc.get_stereo_angle()/180) - A.x()*sin(M_PI*uRwellc.get_stereo_angle()/180));
	XYZ.setZ(A.z());
	
	return XYZ;
}

G4ThreeVector uRwell_strip::intersectionPoint(double m, double c, G4ThreeVector A, G4ThreeVector B){
	
	G4ThreeVector XY;
	double mT = (B.y()-A.y())/(B.x()-A.x());
	double cT = -mT*A.x()+A.y();
	
	XY.setX((cT-c)/(m-mT));
	XY.setY( m*XY.x() +c);
	XY.setZ(A.z());
	


	return XY;
	
}

bool uRwell_strip::pointOnsegment(G4ThreeVector X, G4ThreeVector A, G4ThreeVector B){
	
	
	if((X.x()>= fmin(A.x(), B.x())) && (X.x()<=fmax(A.x(), B.x()))&&(X.y()>=fmin(A.y(), B.y())) && (X.y()<=fmax(A.y(), B.y()))){
		
		return true;
	}else{
		return false;
	}
}

double uRwell_strip::cal_length(G4ThreeVector A, G4ThreeVector B){
	double length=0;
	length = sqrt(pow((A.x()-B.x()),2) + pow((A.y()-B.y()),2)) ;
	return length;
}

double uRwell_strip::GetBinomial(double n, double p){
	double answer;
	answer=CLHEP::RandBinomial::shoot(n,p);
	//Very bad method when n=0 or p close to 0 or 1... return easily -1 in these case.
	//So need to help in the limit condition
	if (answer==-1){
		answer=n;
		if (p==0) answer=0;
	}
	return answer;
}

int uRwell_strip::Number_of_strip(uRwellConstants uRwellc){
	
	int N;
	// C-------------D //
	//  -------------  //
	//   -----------   //
	//    A-------B   //
	
	/*** number of strip in AB***/
	
	int n_AB = abs(2*uRwellc.Xhalf_base/(uRwellc.stripU_pitch/sin(M_PI*uRwellc.get_stereo_angle()/180)));

	/** number of strip in CA **/
	double AC = sqrt((pow((uRwellc.Xhalf_base-uRwellc.Xhalf_Largebase),2) + pow((uRwellc.Yhalf+uRwellc.Yhalf),2)));
    double theta = acos(2*uRwellc.Yhalf/(AC));
	int n_AC = AC/(uRwellc.stripU_pitch/cos(theta-abs(M_PI*uRwellc.get_stereo_angle())/180));

	N = n_AB + n_AC+1;
    return N;
}

int uRwell_strip::strip_id(int i, uRwellConstants uRwell ){
	int ID = 0;
	
	ID = Number_of_strip(uRwell)/2+i;
	
	return ID;
	
}


