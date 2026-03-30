# ODE Tableaus

### Explicit Runge-Kutta Methods

  - `constructEuler` - Euler's 1st order method.
  - `constructHeun()` Heun's order 2 method.
  - `constructRalston()` - Ralston's order 2 method.
  - `constructSSPRK22()` - Explicit SSP method of order 2 using 2 stages.
  - `constructKutta3` - Kutta's classic 3rd order method.
  - `constructSSPRK33()` - Explicit SSP method of order 3 using 3 stages.
  - `constructSSPRK43()` - Explicit SSP method of order 3 using 4 stages.
  - `constructRK4` - The classic 4th order "Runge-Kutta" method.
  - `constructRK438Rule` - The classic 4th order "3/8th's Rule" method.
  - `constructSSPRK104()` - Explicit SSP method of order 4 using 10 stages.
  - `constructBogakiShampine3()` - Bogakai-Shampine's 2/3 method.
  - `constructRKF4()` - Runge-Kutta-Fehlberg 3/4.
  - `constructRKF5()` - Runge-Kutta-Fehlberg 4/5.
  - `constructRungeFirst5()` - Runge's first 5th order method.
  - `constructCassity5()` - Cassity's 5th order method.
  - `constructLawson5()` - Lawson's 5th order method.
  - `constructLutherKonen5` - Luther-Konen's first 5th order method.
  - `constructLutherKonen52()` - Luther-Konen's second 5th order method.
  - `constructLutherKonen53()` - Luther-Konen's third 5th order method.
  - `constructPapakostasPapaGeorgiou5()` - Papakostas and PapaGeorgiou more stable order 5 method.
  - `constructPapakostasPapaGeorgiou52()` - Papakostas and PapaGeorgiou more efficient order 5 method.
  - `constructTsitouras5()` - Tsitouras's order 5 method.
  - `constructBogakiShampine5()` - Bogaki and Shampine's Order 5 method.
  - `constructSharpSmart5()` - Sharp and Smart's Order 5 method.
  - `constructCashKarp()` - Cash-Karp method 4/5.
  - `constructDormandPrince()` - Dormand-Prince 4/5.
  - `constructButcher6()` - Butcher's first order 6 method.
  - `constructButcher62()` - Butcher's second order 6 method.
  - `constructButcher63()` - Butcher's third order 6 method.
  - `constructDormandPrince6()` - Dormand-Prince's 5/6 method.
  - `constructSharpVerner6()` Sharp-Verner's 5/6 method.
  - `constructVerner916()` - Verner's more efficient order 6 method (1991).
  - `constructVerner9162()` - Verner's second more efficient order 6 method (1991).
  - `constructVernerRobust6()` - Verner's "most robust" order 6 method.
  - `constructVernerEfficient6()` - Verner's "most efficient" order 6 method.
  - `constructPapakostas6()` - Papakostas's order 6 method.
  - `constructLawson6()` - Lawson's order 6 method.
  - `constructTsitourasPapakostas6()` - Tsitouras and Papakostas's order 6 method.
  - `constructDormandLockyerMcCorriganPrince6()` - the Dormand-Lockyer-McCorrigan-Prince order 6 method.
  - `constructTanakaKasugaYamashitaYazaki6A()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method A.
  - `constructTanakaKasugaYamashitaYazaki6B()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method B.
  - `constructTanakaKasugaYamashitaYazaki6C()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method C.
  - `constructTanakaKasugaYamashitaYazaki6D()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method D.
  - `constructMikkawyEisa()` - Mikkawy and Eisa's order 6 method.
  - `constructChummund6()` - Chummund's first order 6 method.
  - `constructChummund62()` - Chummund's second order 6 method.
  - `constructHuta6()` - Huta's first order 6 method.
  - `constructHuta62()` - Huta's second order 6 method.
  - `constructVerner6()` - An old order 6 method attributed to Verner.
  - `constructDverk()` - The classic DVERK algorithm attributed to Verner.
  - `constructClassicVerner6()` - A classic Verner order 6 algorithm (1978).
  - `constructButcher7()` - Butcher's order 7 algorithm.
  - `constructClassicVerner7()`- A classic Verner order 7 algorithm (1978).
  - `constructVernerRobust7()` - Verner's "most robust" order 7 algorithm.
  - `constructTanakaYamashitaStable7()` - Tanaka-Yamashita more stable order 7 algorithm.
  - `constructTanakaYamashitaEfficient7()` - Tanaka-Yamashita more efficient order 7 algorithm.
  - `constructSharpSmart7()` - Sharp-Smart's order 7 algorithm.
  - `constructSharpVerner7()` - Sharp-Verner's order 7 algorithm.
  - `constructVerner7()` - Verner's "most efficient" order 7 algorithm.
  - `constructVernerEfficient7()` - Verner's "most efficient" order 7 algorithm.
  - `constructClassicVerner8()` - A classic Verner order 8 algorithm (1978).
  - `constructCooperVerner8()` - Cooper-Verner's first order 8 algorithm.
  - `constructCooperVerner82()` - Cooper-Verner's second order 8 algorithm.
  - `constructTsitourasPapakostas8()` - Tsitouras-Papakostas order 8 algorithm.
  - `constructdverk78()` - The classic order 8 DVERK algorithm.
  - `constructEnrightVerner8()` - Enright-Verner order 8 algorithm.
  - `constructCurtis8()` - Curtis' order 8 algorithm.
  - `constructVerner8()` - Verner's "most efficient" order 8 algorithm.
  - `constructRKF8()` - Runge-Kutta-Fehlberg Order 7/8 method.
  - `constructDormandPrice8()` - Dormand-Prince Order 7/8 method.
  - `constructDormandPrince8_64bit()` - Dormand-Prince Order 7/8 method.
    Coefficients are rational approximations good for 64 bits.
  - `constructVernerRobust9()` - Verner's "most robust" order 9 method.
  - `constructVernerEfficient9()` - Verner's "most efficient" order 9 method.
  - `constructSharp9()` - Sharp's order 9 method.
  - `constructTsitouras9()` - Tsitouras's first order 9 method.
  - `constructTsitouras92()` - Tsitouras's second order 9 method.
  - `constructCurtis10()` - Curtis' order 10 method.
  - `constructOno10()` - Ono's order 10 method.
  - `constructFeagin10Tableau()` - Feagin's order 10 method.
  - `constructCurtis10()` - Curtis' order 10 method.
  - `constructBaker10()` - Baker's order 10 method.
  - `constructHairer10()` Hairer's order 10 method.
  - `constructFeagin12Tableau()` - Feagin's order 12 method.
  - `constructOno12()` - Ono's order 12 method.
  - `constructFeagin14Tableau()` Feagin's order 14 method.

### Implicit Runge-Kutta Methods

  - `constructImplicitEuler` - The 1st order Implicit Euler method.
  - `constructMidpointRule` - The 2nd order Midpoint method.
  - `constructTrapezoidalRule` - The 2nd order Trapezoidal rule (2nd order LobattoIIIA)
  - `constructLobattoIIIA4` - The 4th order LobattoIIIA
  - `constructLobattoIIIB2` - The 2nd order LobattoIIIB
  - `constructLobattoIIIB4` - The 4th order LobattoIIIB
  - `constructLobattoIIIC2` - The 2nd order LobattoIIIC
  - `constructLobattoIIIC4` - The 4th order LobattoIIIC
  - `constructLobattoIIICStar2` - The 2nd order LobattoIIIC*
  - `constructLobattoIIICStar4` - The 4th order LobattoIIIC*
  - `constructLobattoIIID2` - The 2nd order LobattoIIID
  - `constructLobattoIIID4` - The 4th order LobattoIIID
  - `constructRadauIA3` - The 3rd order RadauIA
  - `constructRadauIA5` - The 5th order RadauIA
  - `constructRadauIIA3` - The 3rd order RadauIIA
  - `constructRadauIIA5` - The 5th order RadauIIA

## Tableau Methods

```@docs
DiffEqDevTools.stability_region
OrdinaryDiffEq.ODE_DEFAULT_TABLEAU
```

## Explicit Tableaus

```@docs
DiffEqDevTools.constructEuler
DiffEqDevTools.constructRalston
DiffEqDevTools.constructHeun
DiffEqDevTools.constructKutta3
OrdinaryDiffEq.constructBS3
DiffEqDevTools.constructBogakiShampine3
DiffEqDevTools.constructRK4
DiffEqDevTools.constructRK438Rule
DiffEqDevTools.constructRKF4
DiffEqDevTools.constructRKF5
DiffEqDevTools.constructCashKarp
DiffEqDevTools.constructDormandPrince
OrdinaryDiffEq.constructBS5
DiffEqDevTools.constructPapakostasPapaGeorgiou5
DiffEqDevTools.constructPapakostasPapaGeorgiou52
DiffEqDevTools.constructTsitouras5
DiffEqDevTools.constructLutherKonen5
DiffEqDevTools.constructLutherKonen52
DiffEqDevTools.constructLutherKonen53
DiffEqDevTools.constructRungeFirst5
DiffEqDevTools.constructLawson5
DiffEqDevTools.constructSharpSmart5
DiffEqDevTools.constructBogakiShampine5
DiffEqDevTools.constructCassity5
DiffEqDevTools.constructButcher6
DiffEqDevTools.constructButcher62
DiffEqDevTools.constructButcher63
DiffEqDevTools.constructVernerRobust6
DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6A
DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6B
DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6C
DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6D
DiffEqDevTools.constructHuta6
DiffEqDevTools.constructHuta62
DiffEqDevTools.constructVerner6
DiffEqDevTools.constructDormandPrince6
DiffEqDevTools.constructSharpVerner6
DiffEqDevTools.constructVern6
DiffEqDevTools.constructClassicVerner6
DiffEqDevTools.constructChummund6
DiffEqDevTools.constructChummund62
DiffEqDevTools.constructPapakostas6
DiffEqDevTools.constructLawson6
DiffEqDevTools.constructTsitourasPapakostas6
DiffEqDevTools.constructDormandLockyerMcCorriganPrince6
DiffEqDevTools.constructVernerEfficient6
DiffEqDevTools.constructMikkawyEisa
DiffEqDevTools.constructVernerEfficient7
DiffEqDevTools.constructClassicVerner7
DiffEqDevTools.constructSharpVerner7
DiffEqDevTools.constructTanakaYamashitaStable7
DiffEqDevTools.constructSharpSmart7
DiffEqDevTools.constructTanakaYamashitaEfficient7
DiffEqDevTools.constructVernerRobust7
OrdinaryDiffEq.constructTanYam7
DiffEqDevTools.constructEnrightVerner7
DiffEqDevTools.constructDormandPrince8
DiffEqDevTools.constructRKF8
DiffEqDevTools.constructCooperVerner8
DiffEqDevTools.constructCooperVerner82
DiffEqDevTools.constructTsitourasPapakostas8
DiffEqDevTools.constructEnrightVerner8
DiffEqDevTools.constructdverk78
DiffEqDevTools.constructClassicVerner8
DiffEqDevTools.constructDormandPrince8_64bit
DiffEqDevTools.constructCurtis8
OrdinaryDiffEq.constructTsitPap8
DiffEqDevTools.constructSharp9
DiffEqDevTools.constructTsitouras9
DiffEqDevTools.constructTsitouras92
DiffEqDevTools.constructVernerEfficient9
OrdinaryDiffEq.constructVern9
DiffEqDevTools.constructVerner916
DiffEqDevTools.constructVerner9162
DiffEqDevTools.constructVernerRobust9
DiffEqDevTools.constructFeagin10
DiffEqDevTools.constructFeagin10Tableau
DiffEqDevTools.constructOno10
DiffEqDevTools.constructCurtis10
DiffEqDevTools.constructHairer10
DiffEqDevTools.constructBaker10
DiffEqDevTools.constructFeagin12
DiffEqDevTools.constructOno12
DiffEqDevTools.constructFeagin12Tableau
DiffEqDevTools.constructFeagin14
DiffEqDevTools.constructFeagin14Tableau
```

## Implicit Tableaus

```@docs
DiffEqDevTools.constructImplicitEuler
DiffEqDevTools.constructMidpointRule
DiffEqDevTools.constructTrapezoidalRule
DiffEqDevTools.constructLobattoIIIA4
DiffEqDevTools.constructLobattoIIIB2
DiffEqDevTools.constructLobattoIIIB4
DiffEqDevTools.constructLobattoIIIC2
DiffEqDevTools.constructLobattoIIIC4
DiffEqDevTools.constructLobattoIIICStar2
DiffEqDevTools.constructLobattoIIICStar4
DiffEqDevTools.constructLobattoIIID2
DiffEqDevTools.constructLobattoIIID4
DiffEqDevTools.constructGL2
DiffEqDevTools.constructGL4
DiffEqDevTools.constructGL6
DiffEqDevTools.constructRadauIA3
DiffEqDevTools.constructRadauIA5
DiffEqDevTools.constructRadauIIA3
DiffEqDevTools.constructRadauIIA5
```
