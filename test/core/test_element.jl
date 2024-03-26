@testitem "Element" begin
    @test ElementType(1) == Elements.H
    @test ElementType(2) == Elements.He
    @test ElementType(3) == Elements.Li
    @test ElementType(4) == Elements.Be
    @test ElementType(5) == Elements.B
    @test ElementType(6) == Elements.C
    @test ElementType(7) == Elements.N
    @test ElementType(8) == Elements.O
    @test ElementType(9) == Elements.F
    @test ElementType(10) == Elements.Ne
    @test ElementType(11) == Elements.Na
    @test ElementType(12) == Elements.Mg
    @test ElementType(13) == Elements.Al
    @test ElementType(14) == Elements.Si
    @test ElementType(15) == Elements.P
    @test ElementType(16) == Elements.S
    @test ElementType(17) == Elements.Cl
    @test ElementType(18) == Elements.Ar
    @test ElementType(19) == Elements.K
    @test ElementType(20) == Elements.Ca
    @test ElementType(21) == Elements.Sc
    @test ElementType(22) == Elements.Ti
    @test ElementType(23) == Elements.V
    @test ElementType(24) == Elements.Cr
    @test ElementType(25) == Elements.Mn
    @test ElementType(26) == Elements.Fe
    @test ElementType(27) == Elements.Co
    @test ElementType(28) == Elements.Ni
    @test ElementType(29) == Elements.Cu
    @test ElementType(30) == Elements.Zn
    @test ElementType(31) == Elements.Ga
    @test ElementType(32) == Elements.Ge
    @test ElementType(33) == Elements.As
    @test ElementType(34) == Elements.Se
    @test ElementType(35) == Elements.Br
    @test ElementType(36) == Elements.Kr
    @test ElementType(37) == Elements.Rb
    @test ElementType(38) == Elements.Sr
    @test ElementType(39) == Elements.Y
    @test ElementType(40) == Elements.Zr
    @test ElementType(41) == Elements.Nb
    @test ElementType(42) == Elements.Mo
    @test ElementType(43) == Elements.Tc
    @test ElementType(44) == Elements.Ru
    @test ElementType(45) == Elements.Rh
    @test ElementType(46) == Elements.Pd
    @test ElementType(47) == Elements.Ag
    @test ElementType(48) == Elements.Cd
    @test ElementType(49) == Elements.In
    @test ElementType(50) == Elements.Sn
    @test ElementType(51) == Elements.Sb
    @test ElementType(52) == Elements.Te
    @test ElementType(53) == Elements.I
    @test ElementType(54) == Elements.Xe
    @test ElementType(55) == Elements.Cs
    @test ElementType(56) == Elements.Ba
    @test ElementType(57) == Elements.La
    @test ElementType(58) == Elements.Ce
    @test ElementType(59) == Elements.Pr
    @test ElementType(60) == Elements.Nd
    @test ElementType(61) == Elements.Pm
    @test ElementType(62) == Elements.Sm
    @test ElementType(63) == Elements.Eu
    @test ElementType(64) == Elements.Gd
    @test ElementType(65) == Elements.Tb
    @test ElementType(66) == Elements.Dy
    @test ElementType(67) == Elements.Ho
    @test ElementType(68) == Elements.Er
    @test ElementType(69) == Elements.Tm
    @test ElementType(70) == Elements.Yb
    @test ElementType(71) == Elements.Lu
    @test ElementType(72) == Elements.Hf
    @test ElementType(73) == Elements.Ta
    @test ElementType(74) == Elements.W
    @test ElementType(75) == Elements.Re
    @test ElementType(76) == Elements.Os
    @test ElementType(77) == Elements.Ir
    @test ElementType(78) == Elements.Pt
    @test ElementType(79) == Elements.Au
    @test ElementType(80) == Elements.Hg
    @test ElementType(81) == Elements.Tl
    @test ElementType(82) == Elements.Pb
    @test ElementType(83) == Elements.Bi
    @test ElementType(84) == Elements.Po
    @test ElementType(85) == Elements.At
    @test ElementType(86) == Elements.Rn
    @test ElementType(87) == Elements.Fr
    @test ElementType(88) == Elements.Ra
    @test ElementType(89) == Elements.Ac
    @test ElementType(90) == Elements.Th
    @test ElementType(91) == Elements.Pa
    @test ElementType(92) == Elements.U
    @test ElementType(93) == Elements.Np
    @test ElementType(94) == Elements.Pu
    @test ElementType(95) == Elements.Am
    @test ElementType(96) == Elements.Cm
    @test ElementType(97) == Elements.Bk
    @test ElementType(98) == Elements.Cf
    @test ElementType(99) == Elements.Es
    @test ElementType(100) == Elements.Fm
    @test ElementType(101) == Elements.Md
    @test ElementType(102) == Elements.No
    @test ElementType(103) == Elements.Lr
    @test ElementType(104) == Elements.Rf
    @test ElementType(105) == Elements.Db
    @test ElementType(106) == Elements.Sg
    @test ElementType(107) == Elements.Bh
    @test ElementType(108) == Elements.Hs
    @test ElementType(109) == Elements.Mt
    @test ElementType(110) == Elements.Ds
    @test ElementType(111) == Elements.Rg
    @test ElementType(112) == Elements.Cn
    @test ElementType(113) == Elements.Nh
    @test ElementType(114) == Elements.Fl
    @test ElementType(115) == Elements.Mc
    @test ElementType(116) == Elements.Lv
    @test ElementType(117) == Elements.Ts
    @test ElementType(118) == Elements.Og
    @test ElementType(255) == Elements.Unknown

    @test_throws ArgumentError ElementType(-1)
    @test_throws ArgumentError ElementType(119)

    # see https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/pull/28
    @test_broken parse(ElementType, "MG") == Elements.Mg 
    @test_broken parse(ElementType, "ZN") == Elements.Zn
end
