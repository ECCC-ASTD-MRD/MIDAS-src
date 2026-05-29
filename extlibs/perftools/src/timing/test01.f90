program test
call smp_rebind
stop
end
subroutine smp_rebind
integer i
integer omp_get_thread_num
external omp_get_thread_num
call process_rebind()
!$OMP parallel
print *,'I am thread',omp_get_thread_num()
call thread_rebind()
!$OMP end parallel
return
end
