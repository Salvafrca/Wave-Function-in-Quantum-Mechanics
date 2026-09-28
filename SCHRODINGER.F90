program schrodinger
    implicit none
    real*8 Ax, At, w, L, pi, sigma_0,k_0
    real *8 :: norma , vmx,vmp, vmp_2, vmx_2, Hamiltoniano
    integer  i, j, m,corte
    integer, parameter :: S = 1000
    integer, parameter :: N = 2000
    complex*16, dimension(:),allocatable  :: V, q,fi, fi_,fii_,xdis,fi_aux
    complex*16, dimension(:), allocatable :: alfa, beta , gammia, b
    complex*16 ui,a_2
    ui = (0.d0, 1.d0)
    
    ! Número de intervalos en los que se ve a discretizar el pozo
    m = 3
    ! Longitud del pozo
    L = 1.d0


    allocate(V(0:S), q(0:S), fi(0:S), fi_(0:S), fii_(0:S),xdis(0:S),fi_aux(0:S))
    allocate(alfa(0:S-1), beta(0:S-1), gammia(0:S-1), b(0:S-1))

    ! Por tanto el salto de intervalo es de 
    Ax = L / real(S)
    At = 0.0001d0
    alfa = (0.d0,0.0d0)
    beta = (0.d0,0.0d0)
    pi = 4.d0*atan(1.d0)
    ! Calculo también la w
    w = (At) / ((Ax)*(Ax))
    a_2 = (-2.d0) + (2.d0*ui)/(w)
    corte = 5

    k_0 = 1.d0/(2.d0*sqrt(At))
    sigma_0 = 1.d0/16.d0
    open(10, file='fi_prob.dat', status ='unknown')
    open(11, file='ReIm.dat', status='unknown')
    open(12, file='mediox.dat', status='unknown')
    open(13, file='mediop.dat', status='unknown')
    open(14, file='Hprinc.dat', status='unknown')
    open(15, file='Norm.dat', status='unknown')
    open(16, file='Hamiltoniano.dat', status='unknown')
    

    !do j = 1, S-1
    !    fi(j) = sin(real(m)*pi*j*Ax)
    !end do

    ! FUNCION GAUSSIANA
    do j = 1, S-1
       fi(j) = exp(ui*k_0*j*Ax)* exp(-8.d0*(4*j -real(S))**2)/(real(S)*real(S))
    end do


    norma = sum((abs(fi)**2*Ax))

    do j = 1, S-1
        fi(j) = fi(j)/sqrt(norma)
    end do

    fi(0) = (0.d0,0.d0) 
    fi(S) = (0.d0,0.d0)
    q(0) = (0.d0,0.d0)
    q(S) = (0.d0,0.d0)

    do i = 0, S
        xdis(i) = i*Ax
    end do

    ! La posición se obtendrá viendo cuantos saltos discretos he dado partiendo desde el extremo izquierdo del pozo
    ! y cada salto es una posición, la primera posición es en el origen y se da por tanto tras dar 0 saltos
    ! por tanto j va desde 0 hasta S - 1
    ! Dentro de la caja el potencial es 0
    do i = 0, S
        V(i) = 0.d0
    end do

     !DEFINIR LAS FUNCIONES
     !NORMALIZAR LAS FUNCIOENS DE ONDA INICIAL

    !ALFA
    ! Resolver la ec, eso significa calcular el que???
    ! Evolución de la probabilidad o de la densidad de probabilidad???

    gammia(S-1) = (1.d0)/(a_2)
    alfa(S-1) = (0.d0,0.0d0)

    do j = S-2, 0, -1
        alfa(j) = -gammia(j+1)
        gammia(j) = (1.d0) / (a_2 + alfa(j))
    end do
        



    !BUCLE TEMPORAL
    ! 1º CALCULO BETA CON RECURRENCIA
    ! Cálculo del coeficiente beta para la longitud (j-1) y para el instante n
    ! Primero necesito la gamma
    do j = 1, N
        beta(S-1) = (0.d0,0.0d0)
        do i = S-1 , 0, -1
            b(i) = (4.d0*ui*fi(i)) / (w)
        end do
        
        do i = S-2, 0, -1
            beta(i) = gammia(i+1) * (b(i+1) - beta(i+1))
        end do



        do i = 1, S-1
        ! 2º CALCULO q
        ! Cálculo de q(j + 1,n), la q de un paso necesita la alfa, q, beta del paso anterior
            q(i) = alfa(i-1) * q(i-1) + beta(i-1)
        end do 
        do i = 1, S-1
            ! 3º CALCULO DE fi
            fi_aux(i) = q(i) - fi(i)
        end do

        fi = fi_aux
        fi(0) = (0.d0,0.d0)
        fi(S) = (0.d0,0.d0)
        if (corte == 5) then 
            do i = 0, S
                write(10,*) i*Ax , (abs(fi(i)))**2
                write(11,*) i*Ax, real(fi(i),8), aimag(fi(i))
            end do
            norma = sum((abs(fi)**2*Ax))
            vmx = sum(xdis*abs(fi)**2)*Ax
            vmx_2 = sum(xdis*xdis*abs(fi)**2)*Ax
            do i =1, S-1
                fi_(i)= (fi(i+1)-fi(i))/Ax
                fii_(i)=(fi(i+1)+fi(i-1)-2*fi(i))/(Ax**2)
            end do
        
            vmp = real(sum(conjg(fi(1:S-1))*-(ui)*fi_(1:S-1))*Ax,8)
            vmp_2 = real(sum(conjg(fi(1:S-1))*-(1.d0)*fii_(1:S-1))*Ax,8)
            write(11,*);write(11,*)

            write(10,*);write(10,*)

            write(15,*) j*At, norma 
            write(15,*);write(15,*)

            write(12,*) j*At, vmx
            write(12,*);write(12,*) 

            write(13,*) j*At, vmp 
            write(13,*);write(13,*)

            write(16,*) j*At, vmp_2
            write(16,*);write(16,*)

            write(14,*) j*At, sqrt(vmx_2-vmx*vmx)*sqrt(vmp_2-vmp*vmp)
            write(14,*);write(14,*)
            corte = 0
        end if
        corte = corte +1
    end do
        close(10);close(11);close(12);close(13);close(14);close(15);close(16);
    deallocate(V,q,fi,fi_,fii_,xdis,fi_aux)
    deallocate(alfa,beta,gammia,b)


    




    !------------------------------------------!
    ! CALCULO DE LAS VARIABLES QUE SE ME PIDEN !
    !------------------------------------------!

    ! 1º PARTE REAL DE fi, PARTE IMAGINARIA DE fi PARA CADA INSTANTE POR EL QUE PASO
   

    ! 2º EVOLUCION DE LA PROBABILIDAD DE ENCONTRAR LA PARTICULA A LO LARGO DEL TIEMPO EN FUNCION DE x
    ! SE MUEVE T PERO X QUEDA COMO PARÁMETRO, TENDRE QUE CALCULAR PARA CADA INSTANTE DE TIEMPO 
    ! LA PROBABILIDAD EN CADA X Y LUEGO AVANZO AL SIGUIENTE INSTANTE DE TIEMPO


    ! 3º NORMA DE LA PROBABILIDAD A LO LARGO DEL TIEMPO
    ! PARA CADA INSTANTE DE t calculo el valor de la norma
    ! DEBE DAR EN CADA PASO IGUAL A 1


    ! 4.1º VALOR MEDIO DE LA POSICIÓN, MOMENTO, ENERGIA (SE CONSERVARA) FRENTE AL TIEMPO PARA CADA INSTANTE DE TIEMPO 
    ! 4.2º VALOR MEDIO TEORICO DE LA POSICIÓN, MOMENTO, ENERGIA (SE CONSERVARA) FRENTE AL TIEMPO PARA CADA INSTANTE DE TIEMPO 



    ! 5º CALCULAR EL PRODUCTO AxAp para cada instante de tiempo y compararlo con el teorico 



    

    

    ! 4º Siguiente paso avanzando el indice
    
  
end program schrodinger