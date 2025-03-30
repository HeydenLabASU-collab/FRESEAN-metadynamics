program AVRG
integer::ounit
real(kind=8),dimension(5000):: x1,y1,z1,w1
real(kind=8):: x,y,z,w
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

write(*,'(1x,a)') 'Enter first and last atom number '
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
x1(i)=0.0
y1(i)=0.0
z1(i)=0.0
w1(i)=0.0
enddo

do inp=1,ndata
write(*,'(1x,a,i5,a)')'Enter name of data file',inp
write(*,*)datafile(inp)
open(inunit,file=datafile(inp),status='old')
do i=1,nres
read(inunit,*)x,y,z,w
x1(i)=x1(i)+x
y1(i)=y1(i)+y
z1(i)=z1(i)+z
w1(i)=w1(i)+w
enddo
enddo


do i=1,nres
x1(i)=x1(i)/w1(i)
y1(i)=y1(i)/w1(i)
z1(i)=z1(i)/w1(i)
write(ounit,222)x1(i),y1(i),z1(i)
enddo
222 format (3f8.3)
end

