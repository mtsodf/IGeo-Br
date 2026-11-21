
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE MeshAuxiliares

  CONTAINS

    FUNCTION IsInside(CoordX,CoordY,CoordZ,Nodes,x1,x2,y1,y2,z1,z2)

      IMPLICIT NONE
      LOGICAL:: IsInside
      REAL(8),INTENT(IN):: CoordX(:),CoordY(:),CoordZ(:)
      INTEGER,INTENT(IN):: Nodes(8)
      REAL(8),INTENT(IN):: x1,x2,y1,y2,z1,z2
      
      REAL(8):: x,y,z      
      INTEGER:: num(8),nod=8
      
      IsInside = .FALSE.
      
      x = SUM(CoordX(Nodes))/nod
      y = SUM(CoordY(Nodes))/nod
      z = SUM(CoordZ(Nodes))/nod
      
      IF(z>=z1)THEN
        IF(z<=Z2)THEN
          IF(y>=y1)THEN
            IF(y<=y2)THEN
              IF(x>=x1)THEN
                IF(x<=x2)THEN
                  IsInside = .TRUE.
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
      
    END FUNCTION

    ! Rotina para geração da malha de elementos finitos para um cubo.
    ! Elementos Hexaédricos de 8, 14 e 20 nós.
    ! Smith & Griffiths 2004.
    SUBROUTINE hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)

        ! This subroutine generates nodal coordinates and numbering for
        ! 8, 14 or 20-node "bricks" counting x-z planes in the y-direction.

        IMPLICIT NONE
        INTEGER,INTENT(IN)::iel
        REAL(8),INTENT(IN)::x_coords(:),y_coords(:),z_coords(:)
        REAL(8),INTENT(OUT)::coord(:,:)
        INTEGER,INTENT(OUT)::num(:)
        REAL(8)::pt5=0.5
        INTEGER::fac1,fac2,ip,iq,is,iplane,nod,nxe,nze
        nxe=UBOUND(x_coords,1)-1
        nze=UBOUND(z_coords,1)-1
        nod=UBOUND(num,1)
        iq=(iel-1)/(nxe*nze)+1
        iplane=iel-(iq-1)*nxe*nze
        is=(iplane-1)/nxe+1
        ip=iplane-(is-1)*nxe
        SELECT CASE(nod)
            CASE(8)
                fac1=(nxe+1)*(nze+1)*(iq-1)
                num(1)=fac1+is*(nxe+1)+ip
                num(2)=fac1+(is-1)*(nxe+1)+ip
                num(3)=num(2)+1
                num(4)=num(1)+1
                num(5)=(nxe+1)*(nze+1)*iq+is*(nxe+1)+ip
                num(6)=(nxe+1)*(nze+1)*iq+(is-1)*(nxe+1)+ip
                num(7)=num(6)+1
                num(8)=num(5)+1

                coord(1:2,1)=x_coords(ip)
                coord(5:6,1)=x_coords(ip)
                coord(3:4,1)=x_coords(ip+1)
                coord(7:8,1)=x_coords(ip+1)

                coord(1:4,2)=y_coords(iq)
                coord(5:8,2)=y_coords(iq+1)

                coord(2:3,3)=z_coords(is)
                coord(6:7,3)=z_coords(is)
                coord(1:4:3,3)=z_coords(is+1)
                coord(5:8:3,3)=z_coords(is+1)

            CASE(14)
                fac1=(2*nxe+1)*(2*nze+1)*(iq-1)
                fac2=(2*nxe+1)*(2*nze+1)*iq
                num(1)=fac1+is*(2*nxe+1)+ip
                num(2)=num(1)-(2*nxe+1)
                num(3)=num(2)+1
                num(4)=num(1)+1
                num(5)=num(2)+nxe+1
                num(6)=fac1+(nxe+1)*(nze+1)+nxe*nze+(is-1)*(2*nxe+1)+nxe+ip
                num(7)=num(6)-nxe
                num(8)=num(6)+1
                num(9)=num(8)+nxe
                num(10)=fac2+is*(2*nxe+1)+ip
                num(11)=num(10)-(2*nxe+1)
                num(12)=num(11)+1
                num(13)=num(10)+1
                num(14)=num(11)+nxe+1

                coord(1:2,1)=x_coords(ip)
                coord(6,1)=x_coords(ip)
                coord(10:11,1)=x_coords(ip)
                coord(5:9:2,1)=pt5*(x_coords(ip)+x_coords(ip+1))
                coord(14,1)=pt5*(x_coords(ip)+x_coords(ip+1))
                coord(3:4,1)=x_coords(ip+1)
                coord(8,1)=x_coords(ip+1)
                coord(12:13,1)=x_coords(ip+1)

                coord(1:5,2)=y_coords(iq)
                coord(6:9,2)=pt5*(y_coords(iq)+y_coords(iq+1))
                coord(10:14,2)=y_coords(iq+1)

                coord(2:3,3)=z_coords(is)
                coord(7,3)=z_coords(is)
                coord(11:12,3)=z_coords(is)
                coord(5:6,3)=pt5*(z_coords(is)+z_coords(is+1))
                coord(8:14:6,3)=pt5*(z_coords(is)+z_coords(is+1))
                coord(1:4:3,3)=z_coords(is+1)
                coord(9,3)=z_coords(is+1)
                coord(10:13:3,3)=z_coords(is+1)

            CASE(20)
                fac1=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*(iq-1)
                fac2=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*iq
                num(1)=fac1+(3*nxe+2)*is+2*ip-1
                num(2)=fac1+(3*nxe+2)*is-nxe+ip-1
                num(3)=num(1)-3*nxe-2
                num(4)=num(3)+1
                num(5)=num(4)+1
                num(6)=num(2)+1
                num(7)=num(1)+2
                num(8)=num(1)+1
                num(9)=fac2-(nxe+1)*(nze+1)+(nxe+1)*is+ip
                num(10)=num(9)-nxe-1
                num(11)=num(10)+1
                num(12)=num(9)+1
                num(13)=fac2+(3*nxe+2)*is+2*ip-1
                num(14)=fac2+(3*nxe+2)*is-nxe+ip-1
                num(15)=num(13)-3*nxe-2
                num(16)=num(15)+1
                num(17)=num(16)+1
                num(18)=num(14)+1
                num(19)=num(13)+2
                num(20)=num(13)+1

                coord(1:3,1)=x_coords(ip)
                coord(9:10,1)=x_coords(ip)
                coord(13:15,1)=x_coords(ip)
                coord(5:7,1)=x_coords(ip+1)
                coord(11:12,1)=x_coords(ip+1)
                coord(17:19,1)=x_coords(ip+1)
                coord(4:8:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))
                coord(16:20:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))

                coord(1:8,2)=y_coords(iq)
                coord(13:20,2)=y_coords(iq+1)
                coord(9:12,2)=pt5*(y_coords(iq)+y_coords(iq+1))

                coord(1,3)=z_coords(is+1)
                coord(7:9,3)=z_coords(is+1)
                coord(12:13,3)=z_coords(is+1)
                coord(19:20,3)=z_coords(is+1)
                coord(3:5,3)=z_coords(is)
                coord(10:11,3)=z_coords(is)
                coord(15:17,3)=z_coords(is)
                coord(2:6:4,3)=pt5*(z_coords(is)+z_coords(is+1))
                coord(14:18:4,3)=pt5*(z_coords(is)+z_coords(is+1))

            CASE DEFAULT
                WRITE(11,'(a)')"Wrong number of nodes for hexahedral element"
                STOP
        END SELECT

        RETURN
    END SUBROUTINE hexahedron_xz

    ! Rotina para definição do tamanho da malha de elementos finitos para um cubo.
    ! Smith & Griffiths 2004.
    SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
        !
        !  This subroutine returns the number of elements (nels) and the number
        !  of nodes (nn) in a 2-d geometry-created mesh.
        !
        IMPLICIT NONE
        CHARACTER(LEN=15),INTENT(IN)::element
        INTEGER,INTENT(IN)::nod,nxe,nye
        INTEGER,INTENT(IN),OPTIONAL::nze
        INTEGER,INTENT(OUT)::nels,nn
        IF(element=="triangle")THEN
            nels=nxe*nye*2
            IF(nod==3)nn=(nxe+1)*(nye+1)
            IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
            IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
            IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
        ELSE IF(element=="quadrilateral")THEN
            nels=nxe*nye
            IF(nod==4)nn=(nxe+1)*(nye+1)
            IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
            IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
            IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
        ELSE IF(element=="hexahedron")THEN
            nels=nxe*nye*nze
            IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
            IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
            IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
                (nxe+1)*(nze+1)*nye
        END IF
        RETURN
    END SUBROUTINE mesh_size
    
    SUBROUTINE BuildHexPillars(NX,NY,Lx,Ly,LayersTops,xPillars,yPillars,zPillars)
    
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NX,NY
      REAL(8),INTENT(IN):: Lx,Ly
      REAL(8),INTENT(IN OUT):: LayersTops(:)
      REAL(8),ALLOCATABLE,INTENT(IN OUT):: xPillars(:), yPillars(:), zPillars(:)
      
      INTEGER:: k,i,j,NZ
      
      NZ = SIZE(LayersTops)
      
      ALLOCATE( xPillars(NX+1), yPillars(NY+1), zPillars(NZ+1) )
      xPillars=0.0
      yPillars=0.0
      zPillars=0.0

      !Profundidade das camadas.
      DO k=1 , (NZ)
        zPillars(k) = LayersTops(k)
      END DO
      IF(NZ>1)THEN
        zPillars( NZ+1 ) =  zPillars( NZ ) + zPillars( NZ ) - zPillars( NZ-1 )
      ELSE
        zPillars( NZ+1 ) =  zPillars( NZ ) - 1.0
      END IF
          
      !Criando falsos pilares em x,y para posterior correção de coordenadas.
      DO i=1 , (NX+1)
        xPillars(i) = (i-1)*LX/NX
      END DO
      DO j=1 , (NY+1)
        yPillars(j) = (j-1)*LY/NY
      END DO
    
    END SUBROUTINE
    
  END MODULE MeshAuxiliares
