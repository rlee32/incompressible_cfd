c Takes single-block p3d grid and converts them to a form more easily read by melgen and sefs2d.
	program main
	implicit none
      
      character(len=200),dimension(:),allocatable::filenames
	character(len=200)::gridfile
      integer::i,j,k,nblocks,m,im,jm,km,blockfile
      integer,dimension(:),allocatable::imax,jmax,kmax
      real*8,dimension(:,:,:,:),allocatable::x,y,z

c Read in p3d files and allocate arrays.
      open(unit=12, form='formatted',file='in/blockNames')
	read(12,*) gridfile
	open(unit=11, form='formatted',file='in/'//trim(gridfile))
	read(11,*) nblocks
	print *, "Blocks read in: ",nblocks
	allocate(filenames(nblocks))
	allocate(imax(nblocks))
	allocate(jmax(nblocks))
	allocate(kmax(nblocks))
	read(11,*) ( imax(m),jmax(m),kmax(m), m=1,nblocks )
	im=maxval(imax)
	jm=maxval(jmax)
	km=maxval(kmax)
	print *, "Maximum dimensions read in: ",im,jm,km
	allocate(x(im,jm,km,nblocks))
	allocate(y(im,jm,km,nblocks))
	allocate(z(im,jm,km,nblocks))
	do m=1,nblocks
		read(11,*)
     +	((( x(i,j,k,m),j=1,jmax(m)),i=1,imax(m)),k=1,kmax(m)),
     +	((( y(i,j,k,m),j=1,jmax(m)),i=1,imax(m)),k=1,kmax(m)),
     +	((( z(i,j,k,m),j=1,jmax(m)),i=1,imax(m)),k=1,kmax(m))
		read(12,*) filenames(m)
	enddo
	close(11)

c Write out for elgen (single, multi-block file).
	do m=1,nblocks
	blockfile=20+m
      open ( unit=blockfile, form='formatted',
     +  file = 'out/'//trim(filenames(m))//'.sefs',status='replace')
      write(blockfile,*) imax(m),jmax(m),kmax(m)
      do k=1,kmax(m)
            do j=1,jmax(m)
                  do i=1,imax(m)
      	            write(blockfile,*) x(i,j,k,m),y(i,j,k,m),z(i,j,k,m)
                  enddo
                  write(blockfile,*)
            enddo
            write(blockfile,*)
      enddo
      close(blockfile)
	enddo



	end



c Write out for elgen (single, multi-block file).	
c      open ( unit=12, form='formatted',
c     +  file = 'grid.sefs',status='replace')
c      write(12,*) nblocks
c	do m=1,nblocks
c      write(12,*) imax(m),jmax(m),kmax(m)
c	enddo
c      write(12,*)
c	do m=1,nblocks
c      do k=1,kmax(m)
c            do j=1,jmax(m)
c                  do i=1,imax(m)
c      	            write(12,*) x(i,j,k,m),y(i,j,k,m),z(i,j,k,m)
c                  enddo
c                  write(12,*) 
c            enddo
c            write(12,*) 
c      enddo
c	enddo
c      close(12)
c	end
