
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE MainConstants

    INTEGER::ndim=3,&  !n˙mero de dimensıes do problema
        nod=8,&  !n˙mero de nÛs de deslocamento por elemento
        ndof=24,&  !n˙mero de graus de liberdade de deslocamento por elemento
        nodf=8,&  !n˙mero de nÛs com grau de liberdade de poro-press„o (mas ser· usado como constante).
        nodof=3,&  !n˙mero de graus de liberdade de deslocamento por nÛ
        nst=6  !n˙mero de componentes no tensor de tensıes
    CHARACTER(LEN=320)::element="hexahedron"  !tipo de elemento finito ('hexahedron' ou 'tetrahedron')

  END MODULE
