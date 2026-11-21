
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE MainConstants

    INTEGER::ndim=3,&  !número de dimensões do problema
        nod=8,&  !número de nós de deslocamento por elemento
        ndof=24,&  !número de graus de liberdade de deslocamento por elemento
        nodf=8,&  !número de nós com grau de liberdade de poro-pressão (mas será usado como constante).
        nodof=3,&  !número de graus de liberdade de deslocamento por nó
        nst=6  !número de componentes no tensor de tensões
    CHARACTER(LEN=320)::element="hexahedron"  !tipo de elemento finito ('hexahedron' ou 'tetrahedron')

  END MODULE
