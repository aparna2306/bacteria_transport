module functions
implicit none
contains

function C_in(time) result (res)
double precision, intent(in):: time
double precision:: res

if (time>=0d0 .and. time<(10d0)) then 
res= (1d0*((time-10d0)/(0d0-10d0)))+(0.731d0*(time-0d0)/(10d0-0d0))

elseif (time>=10d0 .and. time<30d0) then
res= (0.731d0*((time-30d0)/(10d0-30d0)))+(0.898d0*(time-10d0)/(30d0-10d0))

elseif (time>=30d0 .and. time<58d0) then
res= (0.898d0*((time-58d0)/(30d0-58d0)))+(1.198d0*(time-30d0)/(58d0-30d0))

elseif (time>=58d0 .and. time<60d0) then
res= (1.198d0*((time-60d0)/(58d0-60d0)))+(1.141d0*(time-58d0)/(60d0-58d0))

elseif (time>=60d0 .and. time<90d0) then
res= (1.141d0*((time-90d0)/(60d0-90d0)))+(0.765d0*(time-60d0)/(90d0-60d0))

elseif (time>=90d0 .and. time<120d0) then
res= (0.765d0*((time-120d0)/(90d0-120d0)))+(0.878d0*(time-90d0)/(120d0-90d0))

elseif (time>=120d0 .and. time<138d0) then
res= (0.878d0*((time-138d0)/(120d0-138d0)))+(1.1368d0*(time-120d0)/(138d0-120d0))

elseif (time>=138d0 .and. time<200d0) then
res= (1.1368d0*((time-200d0)/(138d0-200d0)))+(1.124d0*(time-138d0)/(200d0-138d0))

elseif (time>=200d0 .and. time<240d0) then
res= (1.124d0*((time-240d0)/(200d0-240d0)))+(1.009d0*(time-200d0)/(240d0-200d0))

end if
end function 


end module


