load chimera.pdb
select CaM1, chain A
color cyan, CaM1
center CaM1
label CaM1 and name CA, 'CaM1'
select BLA, chain B
color green, BLA
center BLA
label BLA and name CA, 'BLA'
show cartoon
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.2
show surface
set transparency, 0.7
bg_color white
zoom
