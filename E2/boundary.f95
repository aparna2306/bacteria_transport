module functions
implicit none
double precision, parameter:: pi= 3.1415d0
contains

function v_linear(Q, area, por) result (res)
double precision, intent (in):: Q, area, por
double precision:: res
res= Q/(area*por)
end function

function disp_coeff(dispersion, v_linear) result (res)
double precision, intent (in):: dispersion, v_linear 
double precision:: res
res= dispersion*v_linear
end function

function Peclet(L, dispersion) result (res)
double precision, intent(in):: dispersion, L
double precision:: res
res= L/dispersion 
end function

function pore_volume(L, v_linear) result(res)
double precision, intent(in):: L, v_linear
double precision:: res
res= L/v_linear
end function

function Ret_fac(rho, por, K_a) result(res)
double precision, intent (in):: rho, por, K_a
double precision:: res
res= 1d0 + (rho*K_a/por)
end function

function Da(L, v_linear, K_D) result(res)
double precision, intent (in):: L, v_linear, K_D
double precision:: res
res= L/(v_linear*K_D)
end function

function area(rad) result (res)
double precision, intent(in):: rad
double precision:: res
res= pi*(rad**2d0)
end function


function dNum(u,dx) result(res)
double precision, intent(in):: u,dx
double precision:: res
res= u*dx/2d0
end function
end module
