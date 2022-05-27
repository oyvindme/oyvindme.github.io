// To solve for the structural amplitudes with the wake-oscillator VIV model,
// the package complex.js is used (Copyright (c) 2018 Robert Eisele under MIT license). 
// This package is found at https://github.com/infusion/Complex.js/
// The package is not importen in this file, but in the html file with the plots.
//
// The rest of the code is copyrighted (c) by Ø. M. Ellingsen (2022).
// For usage of the code, contact me on my email
// The above copyright notice and any permission notice shall be included in all
// copies or substantial portions of new software.

function get_Re_and_St(fd, fn, nu) {
  St = 0.2;
  Vcr = fd*fn/St;    // Test Critical velocity
  Re = fd * Vcr / nu; // Test Critical Reynolds number
  // if (Re > 5e5){
  //     St=0.22;
  //     Vcr = fd*fn/St; // New Critical velocity
  //     Re = fd * Vcr / nu; // New Critical Reynolds number
  // }
  return [Re, St];
}

function get_Cd_and_Cl(Re, fh, fd) {
  let Cl = lift_coeff(Re)
  let Cd = drag_coeff(Re, fh/fd)*5/3
  return [Cl, Cd];
}

function drag_coeff(Re, lambda) {
  if (lambda<=6){
    k_a =  0.6+0.129*Math.log10(lambda)
  } else if (lambda<=20){
    k_a = 0.7+0.574*(Math.log10(lambda)-0.778)
  } else{
    k_a = 1
  }
  if (Re<2.5*10**5) {
    Cd0 = 1.2
  } else if (Re <=3.5*10**5) {
    Cd0 = 1.2-3.42*(Math.log10(Re)-5.4)
  } else{
    Cd0 = 0.7
  }
  return Cd0*k_a;
}

function lift_coeff(Re){
    if (Re<2e5){
        Cl = 0.7
    } else if (Re>5e5){
        Cl = 0.2
    } else{
        Cl = 0.7 + (Math.log(Re)-Math.log(2e5))*(0.2-0.7)/(Math.log(5e5)-Math.log(2e5))
    }
    return Cl
}

function VIV_V_and_B (fn, fm, fd, fh, fzeta, fIv, fCl, St, wqs, Re){
  // Function for calculating frequency response using the simplified
  // method of Vickery and Basu.
  // Standard values for St, terrain, rho, nu, Cmode, l and Iv is given
  let aL = 0.4;
  let rho = 1.225; // Density of fluid (standard is air)
  let Cmode = 1/5; // Mode factor
  let l = 1; // Correlation length

  fzeta = fzeta/100
  let Iv = Math.min(0.25, fIv); // Imposed maximum turbulence intensity
  let B  = Iv + 0.1;      // Turbulence bandwidth
  let Vcr = fd*fn/St; // Critical velocity
  let fSc = 2*Math.PI/Math.sqrt(1/fzeta**2-1) * 2 * fm / (rho*fd**2); // Scruton number

  let phi = 1/Math.sqrt(B) * wqs**(3/2) * Math.exp(-0.5 * ((1 - wqs**(-1)) /B)**2);
  let mu  = rho*fd**2/fm; // A mass ratio
  let ga  = Math.sqrt(Math.PI*l)/(2*fh/fd); // Correlation length factor
  let Cf  = fCl /(8 * (Math.PI)**2 * St**2); // A scaled forcing factor

  var Kamax
  var Ka

  // Reynolds determined coefficients from Eurocode
  if (Re <= 10**5){
    Kamax   = 2.0; // Maximum aerodynamic damping coefficient
  } else if (Re>=10**6) {
    Kamax   = 1.0;
  } else if (Re<=5*10**5){  // Interpolation range 10**5 < Re < 5*10**5
    Kamax   = 2.0 + (0.5-2.0) / Math.log10(5) * Math.log10(Re/(10**5));
  } else { // Interpolation range 5*10**5 < Re < 10**6
    Kamax   = 0.5 + (1.0-0.5) / Math.log10(2) * Math.log10(Re/(5*10**5));
  }

  // """Speed test to check if turbulence should be accounted for"""
  // Kvtest = {'0' :  10., 'II' : 10., 'IIIa' : 7., 'IIIb' : 7., 'IV' : 7.}
  // if (Vcr <= Kvtest[fterr]){
  //     Ka = Kamax // Aerodynamic damping with turbulence effect from Eurocode
  //  } else {
  //     Ka = Kamax*(1-3*Iv)
  // }
  if (Vcr <= 7.){
      Ka = Kamax; // Aerodynamic damping with turbulence effect from Eurocode
   } else {
      Ka = Kamax*(1-3*Iv);
  }
  // Solve for amplitude
  let Csh = aL**2 / mu / Ka;
  let b   = fzeta * Csh - aL**2;
  let c   = (Csh * Cf**2 * mu**2 * ga * phi**2)/Cmode;
  let sigma = Math.sqrt(0.5* ( -b + Math.sqrt( b**2 + 4*c ) ) );
  let kp = Math.sqrt(2)*( 1 + 1.2*Math.atan(0.75*(fSc/(4*Math.PI*Ka))**4) );

  return (kp*sigma)
}

function theta_approx(fwq, fM, fD, fA, my_tresh = 10**(-10)){
  // Find x=0 where x=sin^2(theta) and theta is phase difference
  let alpha = fA*fM * fwq**2;
  let beta  = fwq**2 * (1 - fA*fM) - 1;
  let a     = alpha**2*(1 + (fD)**(-2));
  let b     = 2*alpha*beta + 2*alpha - alpha**2/(fD)**2;
  let c     = beta**2 + (fD)**2 - 2*alpha;
  let d     = -1*(fD)**2;
  let delta0= b**2 - 3*a*c;
  let delta1= 2*b**3 - 9*a*b*c + 27*d*a**2;
  let Csqrt = Complex(delta1**2 - 4*delta0**3, 0).sqrt();
  let C     = ((Csqrt.add(delta1)).div(2)).pow(1/3);
  let Xi    = ((Complex(0,1).mul(Math.sqrt(3))).sub(1)).div(2);

  // Solve for x=sin^2(theta), then remove solutions with absolute imaginary part greater than threshhold
  var x1    = ((((Xi.pow(1)).mul(C)).add(b).add(((Xi.pow(-1)).div(C)).mul(delta0))).div(-1*a)).div(3);
  var x2    = ((((Xi.pow(2)).mul(C)).add(b).add(((Xi.pow(-2)).div(C)).mul(delta0))).div(-1*a)).div(3);
  var x3    = ((((Xi.pow(3)).mul(C)).add(b).add(((Xi.pow(-3)).div(C)).mul(delta0))).div(-1*a)).div(3);

  if (Math.abs(x1.im)>my_tresh){
    x1 = NaN;
  } else  {
    x1 = x1.re;
  }
  if (Math.abs(x2.im)>my_tresh){
    x2 = NaN;
  } else  {
    x2 = x2.re;
  }
  if (Math.abs(x3.im)>my_tresh){
    x3 = NaN;
  } else  {
    x3 = x3.re;
  }
  // Find theta values
  var theta1 = Math.asin(x1**(1/2));
  var theta2 = Math.asin(x2**(1/2));
  var theta3 = Math.asin(x3**(1/2));
  if (fwq>1){
    theta1 = Math.PI - theta1;
    theta2 = Math.PI - theta2;
    theta3 = Math.PI - theta3;
  }
  return [theta1, theta2, theta3]
}

function q_approx(fwq, ftheta, fM, fD, fA, feps){
  // evaluate r_q where r_q is amplitude of fluid force multiplied
  my_q = (2*Math.sqrt(1+fA*fM*fwq*Math.sin(ftheta)/(feps*fD)*Math.sin(ftheta)-fD*Math.cos(fwq)));
  return my_q
}

function y_approx(fwq, ftheta, fq, fM, fD){
  // evaluate r_y where r_y is structural amplitude
  my_y = fwq**2*fM/fD*fq*Math.sin(ftheta);
  return my_y
}

function facc_amp(fwq, fM, fD, fA, feps){
  // Find r_y (structural amplitude) without intermediate steps
  // Declare variable
  var theta_1;
  var theta_2;
  var theta_3;
  var my_q1;
  var my_q2;
  var my_q3;
  var my_y1;
  var my_y2;
  var my_y3;
  // Get phase difference
  [theta_1, theta_2, theta_3] = theta_approx(fwq, fM, fD, fA);
  // Get wake amplitude

  if (isNaN(theta_1)){
    my_q1 = NaN;
    my_y1 = NaN;
  } else {
    my_q1 = q_approx(fwq, theta_1, fM, fD, fA, feps);
    my_y1 = y_approx(fwq, theta_1, my_q1, fM, fD);
  }
  if (isNaN(theta_2)){
    my_q2 = NaN;
    my_y2 = NaN;
  } else {
    my_q2 = q_approx(fwq, theta_2, fM, fD, fA, feps);
    my_y2 = y_approx(fwq, theta_2, my_q2, fM, fD);
  }
  if (isNaN(theta_3)){
    my_q3 = NaN;
    my_y3 = NaN;
  } else {
    my_q3 = q_approx(fwq, theta_3, fM, fD, fA, feps);
    my_y3 = y_approx(fwq, theta_3, my_q3, fM, fD);
  }
  return [my_y1, my_y2, my_y3]
}
