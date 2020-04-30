module functions
implicit none
contains

function C_in(time) result (res)
double precision, intent(in):: time
double precision:: res

if (time>=0d0 .and. time<(5d0)) then 
res= (1d0*(time-5d0)/(0d0-5d0))+(0.183d0*(time-0d0)/(5d0-0d0))

else if (time>=5d0 .and. time<(10d0)) then 
res= (0.182d0*(time-10d0)/(5d0-10d0))+(0.142d0*(time-5d0)/(10d0-5d0))

else if (time>=10d0 .and. time<(25d0)) then 
res= (0.142d0*(time-25d0)/(10d0-25d0))+(0.249d0*(time-10d0)/(25d0-10d0))

else if (time>=25d0 .and. time<30d0) then 
res= (0.249d0*(time-30d0)/(25d0-30d0))+(0.105d0*(time-25d0)/(30d0-25d0))

else if (time>=30d0 .and. time<60d0) then 
res= (0.105d0*(time-60d0)/(30d0-60d0))+(0.132d0*(time-30d0)/(60d0-30d0))

else if (time>=60d0 .and. time<70d0) then 
res= (0.132d0*(time-70d0)/(60d0-70d0))+(0.304d0*(time-60d0)/(70d0-60d0))

else if (time>=70d0 .and. time<120d0) then 
res= (0.304d0*(time-120d0)/(70d0-120d0))+(0.101d0*(time-70d0)/(120d0-70d0))

else if (time>=120d0 .and. time<240d0) then 
res= (0.101d0*(time-240d0)/(120d0-240d0))+(0.188d0*(time-120d0)/(240d0-120d0))


end if
end function

end module
