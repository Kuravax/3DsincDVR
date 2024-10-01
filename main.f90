program sinc3d
    implicit none
    integer xnodes, ynodes, znodes, iostatus, memstatus, eigen,io,nodes,i
    double precision, allocatable:: xmatrix(:),ymatrix(:),zmatrix(:), vmatrix(:), hmatrix(:,:), tmatrix(:,:),&
    kmatrix(:),eigenvalues(:), pointmatrix(:,:), umatrix(:,:)
    double precision rminx,rmaxx,rminy,rmaxy,rminz,rmaxz, mass, r0xx,r0yy,r0zz,Fxy,Fxz,Fyz,&
    Fxx,Fyy,Fzz, Lx,Ly,Lz

    open(newunit=io, file="input.txt", status="old", action="read", iostat = iostatus)
    if(iostatus /=0) then
        stop "***File open error***"
    end if

    read(io,*)
    read(io,*) mass, xnodes, ynodes, znodes
    read(io,*)
    read(io,*) rminx,rmaxx
    read(io,*)
    read(io,*) rminy,rmaxy
    read(io,*)
    read(io,*) rminz,rmaxz
    read(io,*)
    read(io,*) r0xx,r0yy, r0zz
    read(io,*)
    read(io,*) Fxx,Fyy,Fzz, Fxy, Fyz, Fxz
    read(io,*)
    read(io,*) eigen
    close(io)

    nodes = xnodes*ynodes*znodes
    Lx = rmaxx-rminx
    Ly = rmaxy-rminy
    Lz = rmaxz-rminz

    allocate (xmatrix(xnodes),ymatrix(ynodes),zmatrix(znodes),tmatrix(nodes,nodes),hmatrix(nodes,nodes),kmatrix(nodes),&
    eigenvalues(nodes),vmatrix(nodes),pointmatrix(4,(5*xnodes+1)*(5*ynodes+1)*(5*znodes+1)),umatrix(4,nodes),STAT = memstatus)
    if(memstatus /= 0) then
        stop "*** Not enough memory for allocation***"
    end if

    call discretegrid(rminx,rmaxx,xnodes,xmatrix)
    call discretegrid(rminy,rmaxy,ynodes,ymatrix)
    call discretegrid(rminz,rmaxz,znodes,zmatrix)
    print*,"Computing Kinetic Energy"
    call threedimkinetic(Kmatrix,xnodes,ynodes,znodes,mass,Lx,Ly,Lz)
    print*,"Computing Potential Energy"
    call threedimpotential(xnodes,ynodes,znodes,xmatrix,ymatrix,zmatrix,vmatrix,&
    Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,r0xx,r0yy,r0zz)
    print*,"Basis transformation..."
    call threedimtransform(xnodes,ynodes,znodes,tmatrix)
    print*,"Computing Hamiltonian"
    call threedimhamiltonian(hmatrix,xnodes,ynodes,znodes,vmatrix,kmatrix,tmatrix)
    print*,"Diagonalizing Hamiltonian"
    call DIAG2(nodes,nodes,eigenvalues,hmatrix)

    open(newunit=io, file="output.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "Computed first 20 eigenvalues of Hamiltonian, Nodes sd", nodes
    write(io,*) "n, energy, hartree "
    do i=1,min(nodes,20)
        write(io,*) i, eigenvalues(i)
        print*, i, eigenvalues(i)
    end do
    write(io,*) "Hamiltonian eigenvectors"
    do i=1,min(nodes,20)
        write(io,*) "Stationary State No. ", i, "Energy", eigenvalues(i)
        write(io,*) hmatrix(:,i)
    end do

    !Print the scatter plot of nth stationary state.
    close(io)
    print*, "Output saved."
    open(newunit=io, file="plots.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "Wavefunction of state n=",eigen
    call FunctionPlot4D(hmatrix,xnodes,ynodes,znodes,pointmatrix,eigen,rminx,&
    rminy,rmaxx,rmaxy,rminz,rmaxz)
        do i=1,(5*xnodes+1)*(5*ynodes+1)*(5*znodes+1)
            write(io,*) pointmatrix(:,i)
        end do
    print*,"Plot saved"
    close(io)

    print*, "Output saved."
    open(newunit=io, file="Uplots.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "U matrix plot"
    call Uplot4D(Vmatrix,xnodes,ynodes,znodes,Umatrix,rminx,rminy,rmaxx,rmaxy,&
rminz,rmaxz)
        do i=1,nodes
            write(io,*) Umatrix(:,i)
        end do
    print*,"Plot saved"
    close(io)

    deallocate(xmatrix,ymatrix,zmatrix,kmatrix,vmatrix,tmatrix, hmatrix, &
    eigenvalues, pointmatrix, umatrix,STAT = memstatus)

    if(memstatus /= 0) then
    stop "*** Deallocation fail***"
    end if

end program

!generate a discrete grid with node points
subroutine discretegrid(r_min,r_max,nodes,xmatrix)
    integer:: nodes, i
    double precision:: r_min, r_max, xmatrix(nodes), step
    step = (r_max-r_min)/(nodes+1)
    do i=1,nodes
    xmatrix(i) = r_min + step*i
    end do
return
end subroutine

!generate kinetic energy operator in three dimensions
subroutine threedimkinetic(Kmat,xnodes,ynodes,znodes,mass,Lx,Ly,Lz)
    integer xnodes, ynodes, znodes, i,j,k
    double precision mass, Lx,Ly,Lz, xconst,yconst,zconst
    double precision Kmat(xnodes*ynodes*znodes)

    xconst =  (acos(-1.0d0)/Lx)**2/(2*mass)
    yconst = (acos(-1.0d0)/Ly)**2/(2*mass)
    zconst =  (acos(-1.0d0)/Lz)**2/(2*mass)

    !K_abc = (pi/Lx)^2 *nx^2/2m + (pi/Ly)^2 *ny^2/2m + (pi/Lz)^2 *nz^2/2m

    do i=1,xnodes
        do j=1, ynodes
            do k=1,znodes
                Kmat((i-1)*xnodes*ynodes+(j-1)*ynodes+k) = xconst*i*i+yconst*j*j+zconst*k*k
            end do
        end do
    end do
return
end subroutine

!Generate 3D potential matrix in DVR representation -- diagonal matrix
subroutine threedimpotential(xnodes,ynodes,znodes, xmatrix,ymatrix,zmatrix,vmatrix,&
    Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,r0xx,r0yy,r0zz)
    integer xnodes, ynodes, znodes, i,j,k
    double precision xmatrix(xnodes),ymatrix(ynodes),zmatrix(znodes),vmatrix(xnodes*ynodes*znodes)
    double precision Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,r0xx,r0yy,r0zz
    do i=1,xnodes
        do j=1,ynodes
            do k=1,znodes
                vmatrix((i-1)*xnodes*ynodes+(j-1)*ynodes+k)= Fxx*(xmatrix(i)-r0xx)**2+Fyy*(ymatrix(j)-r0yy)**2+&
                Fzz*(zmatrix(k)-r0zz)**2
                !vmatrix((i-1)*xnodes*ynodes+(j-1)*ynodes+k) = 3.5d-1
            end do
        end do
    end do
return
end subroutine

!generate 3D transform matrix for sinc basis
subroutine threedimtransform(xnodes,ynodes,znodes,tmatrix)
    integer xnodes, ynodes, znodes, i,j,k,l
    double precision tmatrix(xnodes*ynodes*znodes,xnodes*ynodes*znodes)
    double precision xtemp1,xtemp2,ytemp1,ytemp2,ztemp1,ztemp2,temp,a

    xtemp1 = sqrt(2.0/(xnodes+1))
    xtemp2= acos(-1.0d0)/(xnodes+1)
    ytemp1 = sqrt(2.0/(ynodes+1))
    ytemp2= acos(-1.0d0)/(ynodes+1)
    ztemp1 = sqrt(2.0/(znodes+1))
    ztemp2= acos(-1.0d0)/(znodes+1)
!First compute the matrix of z-transform
    do i=1,znodes
        do j=1,i
            temp = sin(i*j*ztemp2)*ztemp1
            tmatrix(j,i) = temp
            tmatrix(i,j) = temp
        end do
    end do
!Copy z-transform matrix to znodes*ynodes size
    do i=1,ynodes*znodes
        do j=1,znodes*ynodes
            k = mod(i,znodes)
            l = mod(j,znodes)
            if (k==0) then
            k=znodes
            end if
            if (l==0) then
            l=znodes
            end if
            temp = tmatrix(k,l)
            tmatrix(i,j)=temp
        end do
    end do
!Compute direct product with y-transform matrix
    do i=1, ynodes
        do j=1, ynodes
            temp = sin(i*j*ytemp2)*ytemp1
            do k=1, znodes
                do l=1, znodes
                    a = tmatrix((i-1)*znodes+k,(j-1)*znodes+l)
                    tmatrix((i-1)*znodes+k,(j-1)*znodes+l) = a*temp

                end do
            end do
        end do
    end do
!Copy yz-direct product for xnodes*ynodes*znodes size
 do i=1,xnodes*ynodes*znodes
        do j=1,znodes*ynodes*xnodes
            k = mod(i,znodes*ynodes)
            l = mod(j,znodes*ynodes)
            if (k==0) then
            k=znodes*ynodes
            end if
            if (l==0) then
            l=znodes*ynodes
            end if
            temp = tmatrix(k,l)
            tmatrix(i,j)=temp
        end do
    end do
!Compute direct product with x-transform matrix
   do i=1, xnodes
        do j=1, xnodes
            temp = sin(i*j*xtemp2)*xtemp1
            do k=1, znodes*ynodes
                do l=1, znodes*ynodes
                    a = tmatrix((i-1)*znodes*ynodes+k,(j-1)*znodes*ynodes+l)
                    tmatrix((i-1)*znodes*ynodes+k,(j-1)*znodes*ynodes+l) = a*temp
                end do
            end do
        end do
    end do
return
end subroutine

!Generate the three-dimensional hamiltonian matrix
subroutine threedimhamiltonian(hmatrix,xnodes,ynodes,znodes,vmatrix,kmatrix,tmatrix)
    integer xnodes, ynodes,znodes,nodes,i,j,k
    double precision hmatrix(xnodes*ynodes*znodes,xnodes*ynodes*znodes),vmatrix(xnodes*ynodes*znodes),&
    kmatrix(xnodes*ynodes*znodes),tmatrix(xnodes*ynodes*znodes,xnodes*ynodes*znodes)
    double precision temp
    nodes = xnodes*ynodes*znodes
    do i=1,nodes
        Hmatrix(i,i) = Kmatrix(i)
        do j=1, i
            temp=0.0
            do k=1,nodes
                temp = temp + Tmatrix(k,i)*Vmatrix(k)*Tmatrix(j,k)
            end do
            Hmatrix(i,j) = temp + Hmatrix(i,j)
            Hmatrix(j,i) = temp + Hmatrix(j,i)
        end do
    end do
return
end subroutine

subroutine FunctionPlot4D(hmatrix,xnodes,ynodes,znodes,pointmatrix,eigen,rminx,&
rminy,rmaxx,rmaxy,rminz,rmaxz)
    integer i,j,k, xnodes, ynodes,znodes,eigen,l,m,n
    double precision Hmatrix(xnodes*ynodes*znodes,xnodes*ynodes*znodes),&
    pointmatrix(4,(5*xnodes+1)*(5*ynodes+1)*(5*znodes+1))
    double precision rminx,rminy,rmaxx,rmaxy,rminz,rmaxz, stepx, stepy, stepz,&
    const,const1,const2,const3, Lx,Ly,Lz, temp

    Lx = rmaxx-rminx
    Ly = rmaxy-rminy
    Lz = rmaxz-rminz

    stepx = Lx/(5*xnodes)
    stepy = Ly/(5*ynodes)
    stepz = Lz/(5*znodes)
    !phi_ij = 2sqrt(2)/sqrt(LxLyLz) * aij * sin(ipix/Lx)* sin(jpiy/Ly) * sin (kpiz/Lz))
    const = 2.0*sqrt(2.0)/sqrt(Lx*Ly*Lz)
    const1 = acos(-1.0)/Lx
    const2 = acos(-1.0)/Ly
    const3 = acos(-1.0)/Lz
    !Generate coordinate points
    do i=1,(5*xnodes+1)
        do j=1,(5*ynodes+1)
            do k=1,(5*znodes+1)
                pointmatrix(1,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)= rminx+stepx*(i-1)
                pointmatrix(2,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)= rminy+stepy*(j-1)
                pointmatrix(3,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)= rminz+stepz*(k-1)
            end do
        end do
    end do
    !Compute function value
    do i=1,(5*xnodes+1)
        do j=1,(5*ynodes+1)
            do k=1,(5*znodes+1)
                temp = 0.0
                do l=1,xnodes
                    do m=1,ynodes
                        do n=1,znodes
                            temp = temp+const*hmatrix(xnodes*ynodes*(l-1)+ynodes*(m-1)+n,eigen)*&
                    sin(const1*pointmatrix(1,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)*l)*&
                    sin(const2*pointmatrix(2,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)*m)*&
                    sin(const3*pointmatrix(3,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)*n)
                        end do
                    end do
                end do
                pointmatrix(4,(i-1)*(5*ynodes+1)*(5*xnodes+1)+(5*ynodes+1)*(j-1)+k)=temp
            end do
        end do
    end do

return
end subroutine

!Plot the potential energy surface
subroutine Uplot4D(Vmat,xnodes,ynodes,znodes, Upointmatrix,&
    rminx,rminy,rmaxx,rmaxy,rminz,rmaxz)
    integer i,j, xnodes, ynodes, znodes, k
    double precision Vmat(xnodes*ynodes*znodes), Upointmatrix(4,xnodes*ynodes*znodes)
    double precision rminx,rminy,rmaxx,rmaxy,rminz,rmaxz,xstep, ystep,zstep, Lx,Ly,Lz
    Lx = rmaxx-rminx
    Ly = rmaxy-rminy
    Lz = rmaxz-rminz
    xstep = Lx/xnodes
    ystep = Ly/ynodes
    zstep = Lz/znodes
    !Generate coordinate points
    do i=1,xnodes
        do j=1,ynodes
            do k=1, znodes
                Upointmatrix(1,xnodes*ynodes*(i-1)+ynodes*(j-1)+k)= rminx+xstep*(i-1)
                Upointmatrix(2,xnodes*ynodes*(i-1)+ynodes*(j-1)+k)= rminy+ystep*(j-1)
                Upointmatrix(3,xnodes*ynodes*(i-1)+ynodes*(j-1)+k)= rminz+zstep*(k-1)
            end do
        end do
    end do
    !Compute function value
    do i=1,xnodes
        do j=1,ynodes
            do k=1,znodes
                Upointmatrix(4,xnodes*ynodes*(i-1)+ynodes*(j-1)+k) = Vmat(xnodes*ynodes*(i-1)+ynodes*(j-1)+k)
            end do
        end do
    end do
return
end subroutine
