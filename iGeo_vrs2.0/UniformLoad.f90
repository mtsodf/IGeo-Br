
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE UniformLoad


  CONTAINS


    SUBROUTINE TopLoad(GeomechanicalModel,q_load,Load)
        !Constante porepressure.

        USE FEMlib
        USE NumericalLib
        USE DiscretizedModel
        USE NumericalUtilities

        IMPLICIT NONE

        TYPE(Model),INTENT(IN):: GeomechanicalModel
        REAL(8),INTENT(IN):: q_load
        REAL(8),INTENT(IN OUT)::Load(0:)

        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER:: iel,index1,index2,index3  !contadores.

        INTEGER::ndim,&  !dimensıes do problema
            nn,&  !n˙mero de nÛs
            nod,&  !n˙mero de nÛs por elemento
            ndof     !n˙mero de graus de liberdade de deslocamento por elemento.
        REAL(8)::f(3),& !force vector.
            area,n(3),coordsTRI(3,3),area1,area2,normal(3)

        INTEGER,ALLOCATABLE::num(:),&  !conectividade de um elemento.
            nodeDOFs(:)  !vetor com os graus de liberdade de cada elemento.
        REAL(iwp),ALLOCATABLE::coord(:,:),&  !coordenada dos nÛs de um elemento
            g_coord(:,:)  !coordenada dos nÛs da malha

        nn = GeomechanicalModel%NumNodes
        ndim = GeomechanicalModel%NumDimensions
        nod = GeomechanicalModel%NumNodesPerElement
        ndof = GeomechanicalModel%NumDOFUperElement

        !MEMORY ALLOCATION.
        ALLOCATE(num(nod),nodeDOFs(3),coord(nod,ndim),g_coord(ndim,nn))

        !NODAL COORDENATES.
        g_coord=0.0
        g_coord(1,:) = GeomechanicalModel%X
        g_coord(2,:) = GeomechanicalModel%Y
        g_coord(3,:) = GeomechanicalModel%Z

        iel=0
        DO index3=1,GeomechanicalModel%GridDivision(3)             !y
            DO index2=1,GeomechanicalModel%GridDivision(2)         !k (top)
                DO index1=1,GeomechanicalModel%GridDivision(1)     !x (along fracture)
                    iel = iel + 1

                    IF( index2==1 )THEN !IF( index2==GeomechanicalModel%GridDivision(3) )THEN

                        num = GeomechanicalModel%ElementNodes(:,iel)    !g_num(:,iel)
                        coord = TRANSPOSE(g_coord(:,num))

                        !surface area.
                        coordsTRI(1,:)=coord(2,:)
                        coordsTRI(2,:)=coord(3,:)
                        coordsTRI(3,:)=coord(7,:)
                        CALL Triangle( coordsTRI, area1, normal )
                        coordsTRI(1,:)=coord(2,:)
                        coordsTRI(2,:)=coord(7,:)
                        coordsTRI(3,:)=coord(6,:)
                        CALL Triangle( coordsTRI, area2, normal )
                        area = area1+area2

                        !Adding element nodal load to global array.
                        f = 0.0
                        n = 0.0
                        n(3) = -1.0   !z vertical load.
                        f(1) = q_load * area * n(1)
                        f(2) = q_load * area * n(2)
                        f(3) = q_load * area * n(3)

                        nodeDOFs = GeomechanicalModel%NodeDofs(:,num(2))
                        Load(nodeDOFs) = Load(nodeDOFs) + f/4.0
                        nodeDOFs = GeomechanicalModel%NodeDofs(:,num(3))
                        Load(nodeDOFs) = Load(nodeDOFs) + f/4.0
                        nodeDOFs = GeomechanicalModel%NodeDofs(:,num(6))
                        Load(nodeDOFs) = Load(nodeDOFs) + f/4.0
                        nodeDOFs = GeomechanicalModel%NodeDofs(:,num(7))
                        Load(nodeDOFs) = Load(nodeDOFs) + f/4.0

                    END IF

                END DO
            END DO
        END DO

        Load(0) = 0.0

        RETURN
    END SUBROUTINE

END MODULE
