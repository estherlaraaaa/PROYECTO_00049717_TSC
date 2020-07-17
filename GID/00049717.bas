*GenData(1) *GenData(2) *GenData(3) *GenData(4) *GenData(5) *GenData(6)
*Set Cond Dirichlet_U *nodes
*Set var NDU=CondNumEntities(int)
*Set Cond Dirichlet_V *nodes
*Set var NDV=CondNumEntities(int)
*Set Cond Dirichlet_W *nodes
*Set var NDW=CondNumEntities(int)
*Set Cond Dirichlet_P *nodes
*set var NDP=CondNumEntities(int)
*npoin *nelem *NDU *NDV *NDW *NDP

Coordinates
*set elems(all)
*loop nodes
*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real)
*end nodes
EndCoordinates

Elements
*loop elems
*ElemsNum *ElemsConec
*end elems
EndElements

Dirichlet_U
*Set Cond Dirichlet_U *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(u,real)
*end nodes
EndDirichlet

Dirichlet_V
*Set Cond Dirichlet_V *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(v,real)
*end nodes
EndDirichlet

Dirichlet_W
*Set Cond Dirichlet_W *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(w,real)
*end nodes
EndDirichlet

Dirichlet_P
*Set Cond Dirichlet_P *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(p,real)
*end nodes
EndDirichlet
