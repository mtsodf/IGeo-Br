
    !GEOMECANICA COMPUTACIONAL DE RESERVATORIOS
    !PETROBRAS-E&P-ENGP / RR / Eng. de Reservatórios
    !Eng. M.Sc. Luis Carlos de Sousa Jr
    !cxnx / sousajr@petrobras.com.br
    !
    !versão 2.0 : 01 de stembro de 2010
    !última modificação: 11 de janeiro de 2012
    
    !****************************************************************************
    !             ROTINAS AUXILIARES PARA PROCESSAMENTO PARALELO
    !****************************************************************************
    
    MODULE ParallelRoutines
  
    CONTAINS
  
    SUBROUTINE PARTICIONAMENTO(vetor,particoes)
      !Particionamento do domínio de acordo com os valores de um vetor de carga.
      !Regra: soma acumulada dos valores absolutos do vetor.
      !
      !VARIÁVEIS
        !IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,cont,&  !contadores
                 nParticoes,&
                 neq
        REAL(iwp)::soma,acumulado
        INTEGER,INTENT(IN OUT)::particoes(:,:)
        REAL(iwp),INTENT(IN)::vetor(:)
      !
      !CORPO DA ROTINA
      !
        nParticoes = SIZE( particoes(1,:) )
        neq = particoes(2,nParticoes)
        !
        soma = SUM( ABS(vetor) )
        acumulado=0.0_iwp
        cont = 0
        !
        DO i=1,neq
          acumulado = acumulado + ABS( vetor(i) )
          IF( acumulado > (cont+1)*soma/nParticoes ) THEN
            cont = cont + 1
            particoes(2,cont) = i
          END IF
        END DO
        !
        DO i=2,nParticoes
          particoes(1,i) = particoes(2,i-1) + 1
        END DO
        !
        !ajustes
          particoes(1,1) = 1
          particoes(2,nParticoes) = neq
        !
    END SUBROUTINE
    
    FUNCTION DOT_PRODUCT_pp(a,b,particoes)
      !Rotina para fazer produto escalar entre vetores.
      ! a . b = SUM( ai*bi )
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp)::DOT_PRODUCT_pp
        REAL(iwp),ALLOCATABLE::res(:)
        REAL(iwp),INTENT(IN)::a(:),b(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      ALLOCATE( res(nParticoes) )
      !
      res = 0.0_iwp
      !

    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        !
        res(i) = DOT_PRODUCT( a(posicao1:posicao2), b(posicao1:posicao2) )
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
      DOT_PRODUCT_pp = SUM(res)
      !
    RETURN
    END FUNCTION DOT_PRODUCT_pp
    
    SUBROUTINE SpAx(a,ia,ja,i1,i2,x,y)
      !Multiplica uma matriz esparsa no formato CSR por um vetor.
      !y = A * x
      !
      !VARIÁVEIS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,&  !contadores
                 n,&  !tamanho dos vetores que serão multiplicados
                 neq,&  !dimensão da matriz (neq x neq)
                 posicao1,posicao2  !marcadores de posição na matriz "a" (CSR)
        REAL(iwp)::ddoti  !MKL - BLAS: função para fazer o dot product de um vetor esparso por um vetor cheio
        INTEGER,INTENT(IN)::i1,i2
        INTEGER,INTENT(IN)::ia(:),ja(:)
        REAL(iwp),INTENT(IN)::a(:)
        REAL(iwp),INTENT(IN)::x(:)
        REAL(iwp),INTENT(OUT)::y(:)
      !
      !CORPO DA ROTINA
      !
      neq = SIZE(ia) - 1
      !
    !DEC$ IVDEP : LOOP
      DO i=i1,i2
        !
        posicao1 = ia(i)
        posicao2 = ia(i+1) - 1
        n = posicao2 - posicao1 + 1
        !
        y(i) = ddoti(n, a( posicao1:posicao2 ), ja( posicao1:posicao2 ), x)
        !
      END DO
      !
    RETURN
    END SUBROUTINE
    
    SUBROUTINE MAX_pp(a,b,particoes,max,maxDif)
      !Acha máximo valor absoluto de a_pp e a máxima diferença absoluta 
      !entre a_pp e b_pp.
      !
      !VARIÁVEIS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN OUT)::a(:),b(:)
        REAL(iwp),INTENT(OUT)::max,maxDif
        REAL(iwp),ALLOCATABLE::temp1(:),temp2(:)
        REAL(iwp)::temp
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      max=0.0
      maxDif=0.0
      ALLOCATE( temp1(nParticoes), temp2(nParticoes) )
      temp1=0.0
      temp2=0.0
      !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i,posicao1,posicao2,temp)
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        !
        temp =  MAXVAL( ABS( a( posicao1:posicao2 ) ) ) 
        !
        IF(temp>temp1(i)) temp1(i) = temp
        !
        temp =  MAXVAL( ABS( a( posicao1:posicao2 ) - b( posicao1:posicao2 ) ) )
        !
        IF(temp>temp2(i)) temp2(i) = temp
        !
      END DO
!$OMP END PARALLEL DO
      max = MAXVAL(temp1)
      maxDif = MAXVAL(temp2)
      ! 
    RETURN
    END SUBROUTINE MAX_pp
    
    SUBROUTINE ADICAOaaeb_pp(a,b,escalar,particoes)
      !Rotina para somar dois vetores e guardar resultado em um terceiro vetor.
      ! a = a + escalar * b
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2,n
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN)::escalar
        REAL(iwp),INTENT(IN)::b(:)
        REAL(iwp),INTENT(IN OUT)::a(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      !
    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2,n) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        n = posicao2 - posicao1 + 1
        !
        !a = a + escalar * b
          a( posicao1:posicao2 ) = a( posicao1:posicao2 ) + escalar * b( posicao1:posicao2 )
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
    RETURN
    END SUBROUTINE ADICAOaaeb_pp
    
    SUBROUTINE ADICAOabea_pp(a,b,escalar,particoes)
      !Rotina para somar dois vetores e guardar resultado em um terceiro vetor.
      ! a = b + escalar * a
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN)::escalar
        REAL(iwp),INTENT(IN)::b(:)
        REAL(iwp),INTENT(IN OUT)::a(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      !
    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        !
          a( posicao1:posicao2 ) = b( posicao1:posicao2 ) + escalar * a( posicao1:posicao2 )
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
    RETURN
    END SUBROUTINE ADICAOabea_pp
    
    SUBROUTINE MULTab_pp(c,a,b,particoes)
      !Rotina para multiplicar dois vetores.
      ! c = a * b
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN)::a(:),b(:)
        REAL(iwp),INTENT(IN OUT)::c(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      !
    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        !
          c( posicao1:posicao2 ) = a( posicao1:posicao2 ) * b( posicao1:posicao2 )
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
    RETURN
    END SUBROUTINE MULTab_pp
    
    SUBROUTINE MULTaeb_pp(a,b,escalar,particoes)
      !Rotina para somar dois vetores e guardar resultado em um terceiro vetor.
      ! a = escalar * b
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN)::escalar
        REAL(iwp),INTENT(IN)::b(:)
        REAL(iwp),INTENT(IN OUT)::a(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      !
    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        !
          a( posicao1:posicao2 ) = escalar * b( posicao1:posicao2 )
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
    RETURN
    END SUBROUTINE MULTaeb_pp
    
    SUBROUTINE SWAP_pp(a,b,particoes)
      !Rotina para somar dois vetores e guardar resultado em um terceiro vetor.
      ! a = escalar * b
      !
      USE OMP_LIB  !bibliotecas OpenMP(R)
      !VARIÁVELS
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i,nParticoes,posicao1,posicao2,N
        INTEGER,INTENT(IN)::particoes(:,:)
        REAL(iwp),INTENT(IN)::b(:)
        REAL(iwp),INTENT(IN OUT)::a(:)
        REAL(iwp),ALLOCATABLE::c(:)
      !
      !CORPO DA ROTINA
      !
      nParticoes = SIZE( particoes(1,:) )
      ALLOCATE( c( size(a) ) )
      !
    !$OMP PARALLEL
    !$OMP DO PRIVATE(posicao1,posicao2,n) SCHEDULE(STATIC)
    !DEC$ IVDEP : LOOP
      DO i=1,nParticoes
        !
        posicao1 = particoes(1,i)
        posicao2 = particoes(2,i)
        n = posicao2 - posicao1 + 1
        !
        !swap
          CALL dswap(n, a( posicao1:posicao2 ), 1, b( posicao1:posicao2 ), 1) 
        !
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
      !
    RETURN
    END SUBROUTINE SWAP_pp
  
  END MODULE