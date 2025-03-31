program AVRG
integer::ounit
real,dimension(300,300):: dccm
character(len=100),dimension(25)::datafile
character(len=100):: outfile
logical:: there,ANS
data inunit,ounit /11,12/

!---------------------------------------------
write(*,'(1x,a)')'Enter the name of output file'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------

write(*,'(1x,a)') 'Enter first and last res number '
read(*,*)ifr,ilr
write(*,*)ifr,ilr
nres=ilr-ifr+1
write(*,*)'total number of res',nres

!write(*,'(1x,a)')'Enter the number of frames in each file'
!read(*,*) nframes
!write(*,*)'# of frames',nframes

write(*,'(1x,a)')'Enter the number of data file'
read(*,*)ndata

write(*,*)'# of data files',ndata

do inp=1,ndata
write(*,'(1x,a,i5,a)')'Enter name of data file',inp
read(*,'(a)')datafile(inp)
write(*,*)datafile(inp)
inquire(file=datafile(inp),exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: INP FILE NOT FOUND'
stop
endif
enddo
do i=1,nres
do j=1,nres
dccm(i,j)=0.0
enddo
enddo

do inp=1,ndata
write(*,'(1x,a,i5,a)')'Enter name of data file',inp
write(*,*)datafile(inp)
open(inunit,file=datafile(inp),status='old')
do i=1,nres
do j=1,nres
read(inunit,*)irj,irk,dcc
dccm(i,j)=dccm(i,j)+dcc
enddo
read(inunit,*)
enddo
enddo


do i=1,nres
do j=1,nres
dccm(i,j)=dccm(i,j)/ndata
write(ounit,111)i,j,dccm(i,j)
111 format(i4,i4,f10.4)
enddo
write(ounit,*)
enddo
end

