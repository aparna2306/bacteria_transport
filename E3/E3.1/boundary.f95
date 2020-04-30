module functions
implicit none
contains

function C_in(time) result (res)
double precision, intent(in):: time
double precision:: res
if (time>=0d0 .and. time<(5d0)) then 
res= (1d0*(time-5d0)/(0d0-5d0))+(0.058064d0*(time-0d0)/(5d0-0d0))

elseif (time>=5d0 .and. time<10d0) then
res= (0.058064d0*(time-10d0)/(5d0-10d0))+(0.028d0*(time-5d0)/(10d0-5d0))

elseif (time>=10d0 .and. time<15d0) then
res= (0.028d0*(time-15d0)/(10d0-15d0))+(0.021d0*(time-10d0)/(15d0-10d0))

elseif (time>=15d0 .and. time<30d0) then
res= (0.021d0*(time-30d0)/(15d0-30d0))+(0.0089d0*(time-15d0)/(30d0-15d0))

elseif (time>=30d0 .and. time<60d0) then
res= (0.0089d0*(time-60d0)/(30d0-60d0))+(0.012d0*(time-30d0)/(60d0-30d0))

elseif (time>=60d0 .and. time<120d0) then
res= (0.012d0*(time-120d0)/(60d0-120d0))+(0.0082d0*(time-60d0)/(120d0-60d0))

elseif (time>=120d0 .and. time<240d0) then
res= (0.0082d0*(time-240d0)/(120d0-240d0))+(0.0114d0*(time-120d0)/(240d0-120d0))
end if
end function 

end module