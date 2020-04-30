module functions
implicit none
contains

function C_in(time) result (res)
double precision, intent(in):: time
double precision:: res

if (time>=0d0 .and. time<(10d0)) then 
res= (1d0*(time-10d0)/(0d0-10d0))+(0.915d0*(time-0d0)/(10d0-0d0))

else if (time>=10d0 .and. time<(15d0)) then 
res= (0.915d0*(time-15d0)/(10d0-15d0))+(0.864d0*(time-10d0)/(15d0-10d0))

else if (time>=15d0 .and. time<(30d0)) then 
res= (0.864d0*(time-30d0)/(15d0-30d0))+(1.444d0*(time-15d0)/(30d0-15d0))

else if (time>=30d0 .and. time<60d0) then
res= (1.444d0*(time-60d0)/(30d0-60d0))+(1.173d0*(time-30d0)/(60d0-30d0))

else if (time>=60d0 .and. time<120d0) then 
res= (1.173d0*(time-120d0)/(60d0-120d0))+(1.752d0*(time-60d0)/(120d0-60d0))

else if (time>=120d0 .and. time<145d0) then 
res= (1.752d0*(time-145d0)/(120d0-145d0))+(2.687d0*(time-120d0)/(145d0-120d0))

else if (time>=145d0 .and. time<240d0) then 
res= (2.687d0*(time-240d0)/(145d0-240d0))+(4.280d0*(time-145d0)/(240d0-145d0))


end if
end function

end module







