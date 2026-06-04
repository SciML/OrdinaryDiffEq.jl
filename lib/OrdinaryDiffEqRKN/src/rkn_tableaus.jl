abstract type NystromConstantCache <: OrdinaryDiffEqConstantCache end

"""
    NystromVITableau{T, T2}

Tableau for velocity-independent Nyström methods.
Fields:
- `a`: nstages × nstages lower-triangular position coupling matrix
- `b`: position update weights (length nstages; b[end] may be 0 if last stage not in position update)
- `bp`: velocity update weights (length nstages)
- `btilde`: embedded position error weights (empty if non-adaptive)
- `bptilde`: embedded velocity error weights (empty if non-adaptive)
- `c`: time nodes for stages 2..nstages (length nstages-1); c[i] is node for stage i+1
- `pos_only_error`: if true, error estimate uses only position components (ERKN5 behaviour)
"""
struct NystromVITableau{T, T2}
    a::Matrix{T}
    b::Vector{T}
    bp::Vector{T}
    btilde::Vector{T}
    bptilde::Vector{T}
    c::Vector{T2}
    pos_only_error::Bool
    # Length of integrator.k slot used by the dense-output interpolant.
    # 2 for Hermite (standard); methods with custom continuous output (e.g. DPRKN6)
    # use 3 to cache the extra stage derivatives.
    kshortsize::Int
end

function DPRKN4Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 4, 4)
    a[2, 1] = convert(T, 1 // 32)
    a[3, 1] = convert(T, 7 // 1000); a[3, 2] = convert(T, 119 // 500)
    a[4, 1] = convert(T, 1 // 14); a[4, 2] = convert(T, 8 // 27); a[4, 3] = convert(T, 25 // 189)
    b = [convert(T, 1 // 14), convert(T, 8 // 27), convert(T, 25 // 189), zero(T)]
    bp = [convert(T, 1 // 14), convert(T, 32 // 81), convert(T, 250 // 567), convert(T, 5 // 54)]
    btilde = [
        convert(T, 1 // 14 + 7 // 150), convert(T, 8 // 27 - 67 // 150),
        convert(T, 25 // 189 - 3 // 20), convert(T, 1 // 20),
    ]
    bptilde = [
        convert(T, 1 // 14 - 13 // 21), convert(T, 32 // 81 + 20 // 27),
        convert(T, 250 // 567 - 275 // 189), convert(T, 5 // 54 + 1 // 3),
    ]
    c = T2[convert(T2, 1 // 4), convert(T2, 7 // 10), one(T2)]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function DPRKN5Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 6, 6)
    a[2, 1] = convert(T, 1 // 128)
    a[3, 1] = convert(T, 1 // 96); a[3, 2] = convert(T, 1 // 48)
    a[4, 1] = convert(T, 1 // 24); a[4, 3] = convert(T, 1 // 12)
    a[5, 1] = convert(T, 9 // 128); a[5, 3] = convert(T, 9 // 64); a[5, 4] = convert(T, 9 // 128)
    a[6, 1] = convert(T, 7 // 90); a[6, 3] = convert(T, 4 // 15); a[6, 4] = convert(T, 1 // 15)
    a[6, 5] = convert(T, 4 // 45)
    b = [
        convert(T, 7 // 90), zero(T), convert(T, 4 // 15), convert(T, 1 // 15),
        convert(T, 4 // 45), zero(T),
    ]
    bp = [
        convert(T, 7 // 90), zero(T), convert(T, 16 // 45), convert(T, 2 // 15),
        convert(T, 16 // 45), convert(T, 7 // 90),
    ]
    btilde = [
        convert(T, 7 // 90 - 1 // 6), zero(T), convert(T, 4 // 15),
        convert(T, 1 // 15 - 1 // 3), convert(T, 4 // 45), zero(T),
    ]
    bptilde = [
        convert(T, 7 // 90), zero(T), convert(T, 16 // 45 - 2 // 3),
        convert(T, 2 // 15 + 1 // 3), convert(T, 16 // 45 - 2 // 3),
        convert(T, 7 // 90),
    ]
    c = T2[
        convert(T2, 1 // 8), convert(T2, 1 // 4), convert(T2, 1 // 2),
        convert(T2, 3 // 4), one(T2),
    ]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function DPRKN6FMTableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 6, 6)
    a[2, 1] = convert(T, 1 // 200)
    a[3, 1] = convert(T, -1 // 2200); a[3, 2] = convert(T, 1 // 22)
    a[4, 1] = convert(T, 637 // 6600); a[4, 2] = convert(T, -7 // 110); a[4, 3] = convert(T, 7 // 33)
    a[5, 1] = convert(T, 225437 // 1968750); a[5, 2] = convert(T, -30073 // 281250)
    a[5, 3] = convert(T, 65569 // 281250); a[5, 4] = convert(T, -9367 // 984375)
    a[6, 1] = convert(T, 151 // 2142); a[6, 2] = convert(T, 5 // 116)
    a[6, 3] = convert(T, 385 // 1368); a[6, 4] = convert(T, 55 // 168); a[6, 5] = convert(T, -6250 // 28101)
    b = [
        convert(T, 151 // 2142), convert(T, 5 // 116), convert(T, 385 // 1368),
        convert(T, 55 // 168), convert(T, -6250 // 28101), zero(T),
    ]
    bp = [
        convert(T, 151 // 2142), convert(T, 25 // 522), convert(T, 275 // 684),
        convert(T, 275 // 252), convert(T, -78125 // 112404), convert(T, 1 // 12),
    ]
    btilde = [
        convert(T, 151 // 2142 - 1349 // 157500), convert(T, 5 // 116 - 7873 // 50000),
        convert(T, 385 // 1368 - 192199 // 900000), convert(T, 55 // 168 - 521683 // 2100000),
        convert(T, -6250 // 28101 + 16 // 125), zero(T),
    ]
    bptilde = [
        convert(T, 151 // 2142 - 1349 // 157500), convert(T, 25 // 522 - 7873 // 45000),
        convert(T, 275 // 684 - 27457 // 90000), convert(T, 275 // 252 - 521683 // 630000),
        convert(T, -78125 // 112404 + 2 // 5), zero(T),
    ]
    c = T2[
        convert(T2, 1 // 10), convert(T2, 3 // 10), convert(T2, 7 // 10),
        convert(T2, 17 // 25), one(T2),
    ]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function DPRKN8Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 9, 9)
    a[2, 1] = convert(T, 1 // 800)
    a[3, 1] = convert(T, 1 // 600); a[3, 2] = convert(T, 1 // 300)
    a[4, 1] = convert(T, 9 // 200); a[4, 2] = convert(T, -9 // 100); a[4, 3] = convert(T, 9 // 100)
    a[5, 1] = convert(T, -66701 // 197352); a[5, 2] = convert(T, 28325 // 32892)
    a[5, 3] = convert(T, -2665 // 5482); a[5, 4] = convert(T, 2170 // 24669)
    a[6, 1] = convert(T, 22701_5747 // 30425_1000); a[6, 2] = convert(T, -5489_7451 // 30425_100)
    a[6, 3] = convert(T, 12942_349 // 10141_700); a[6, 4] = convert(T, -9499 // 304_251)
    a[6, 5] = convert(T, 539 // 9250)
    a[7, 1] = convert(T, -11318_91597 // 9017_89000); a[7, 2] = convert(T, 4196_4921 // 1288_2700)
    a[7, 3] = convert(T, -6663_147 // 3220_675); a[7, 4] = convert(T, 270_954 // 644_135)
    a[7, 5] = convert(T, -108 // 5875); a[7, 6] = convert(T, 114 // 1645)
    a[8, 1] = convert(T, 138_369_59 // 3667458); a[8, 2] = convert(T, -177_314_50 // 1833729)
    a[8, 3] = convert(T, 106_3919_505 // 15647_8208); a[8, 4] = convert(T, -332_138_45 // 3911_9552)
    a[8, 5] = convert(T, 133_35 // 285_44); a[8, 6] = convert(T, -705 // 14272); a[8, 7] = convert(T, 1645 // 57088)
    a[9, 1] = convert(T, 223 // 7938); a[9, 3] = convert(T, 1175 // 8064); a[9, 4] = convert(T, 925 // 6048)
    a[9, 5] = convert(T, 41 // 448); a[9, 6] = convert(T, 925 // 14112); a[9, 7] = convert(T, 1175 // 72576)
    b = [
        convert(T, 223 // 7938), zero(T), convert(T, 1175 // 8064), convert(T, 925 // 6048),
        convert(T, 41 // 448), convert(T, 925 // 14112), convert(T, 1175 // 72576), zero(T), zero(T),
    ]
    bp = [
        convert(T, 223 // 7938), zero(T), convert(T, 5875 // 36288), convert(T, 4625 // 21168),
        convert(T, 41 // 224), convert(T, 4625 // 21168), convert(T, 5875 // 36288),
        convert(T, 223 // 7938), zero(T),
    ]
    btilde = [
        convert(T, 223 // 7938 - 7987_313 // 10994_1300), zero(T),
        convert(T, 1175 // 8064 - 1610_737 // 4467_4560),
        convert(T, 925 // 6048 - 10023_263 // 3350_5920),
        convert(T, 41 // 448 + 497_221 // 1240_9600),
        convert(T, 925 // 14112 - 1002_3263 // 7818_0480),
        convert(T, 1175 // 72576 - 1610_737 // 40207_1040), zero(T), zero(T),
    ]
    bptilde = [
        convert(T, 223 // 7938 - 7987_313 // 10994_1300), zero(T),
        convert(T, 5875 // 36288 - 1610_737 // 4020_7104),
        convert(T, 4625 // 21168 - 1002_3263 // 2345_4144),
        convert(T, 41 // 224 + 497_221 // 620_4800),
        convert(T, 4625 // 21168 - 1002_3263 // 2345_4144),
        convert(T, 5875 // 36288 - 1610_737 // 40207_104),
        convert(T, 223 // 7938 + 4251_941 // 5497_0650), convert(T, -3 // 20),
    ]
    c = T2[
        convert(T2, 1 // 20), convert(T2, 1 // 10), convert(T2, 3 // 10), convert(T2, 1 // 2),
        convert(T2, 7 // 10), convert(T2, 9 // 10), one(T2), one(T2),
    ]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function DPRKN12Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 17, 17)
    a[2, 1] = convert(T, 1 // 5000)
    a[3, 1] = convert(T, 1 // 3750); a[3, 2] = convert(T, 1 // 1875)
    a[4, 1] = convert(T, 7 // 2400); a[4, 2] = convert(T, -1 // 240); a[4, 3] = convert(T, 1 // 160)
    a[5, 1] = convert(T, 2 // 1215); a[5, 3] = convert(T, 4 // 729); a[5, 4] = convert(T, 32 // 18225)
    a[6, 1] = convert(T, 152 // 78125); a[6, 3] = convert(T, 1408 // 196875)
    a[6, 4] = convert(T, 2048 // 703125); a[6, 5] = convert(T, 432 // 546875)
    a[7, 1] = convert(T, 29 // 51200); a[7, 3] = convert(T, 341 // 387072)
    a[7, 4] = convert(T, -151 // 345600); a[7, 5] = convert(T, 243 // 716800)
    a[7, 6] = convert(T, -11 // 110592)
    a[8, 1] = convert(T, 37 // 12000); a[8, 4] = convert(T, 2 // 1125)
    a[8, 5] = convert(T, 27 // 10000); a[8, 6] = convert(T, 5 // 3168); a[8, 7] = convert(T, 224 // 20625)
    a[9, 1] = convert(T, 100467472123373 // 27511470744477696)
    a[9, 3] = convert(T, 101066550784375 // 25488568483854336)
    a[9, 4] = convert(T, 49478218404275 // 15475202293768704)
    a[9, 5] = convert(T, 21990175014231 // 2674726322379776)
    a[9, 6] = convert(T, -3576386017671875 // 2723635603703291904)
    a[9, 7] = convert(T, 16163228153 // 1654104722787)
    a[9, 8] = convert(T, 38747524076705 // 10316801529179136)
    a[10, 1] = convert(T, 62178936641284701329 // 16772293867250014666848)
    a[10, 3] = convert(T, 46108564356250 // 9072835168325103)
    a[10, 4] = convert(T, 1522561724950 // 1296119309760729)
    a[10, 5] = convert(T, -45978886013453735443 // 2174186242050927827184)
    a[10, 6] = convert(T, 299403512366617849203125 // 4981371278573254356053856)
    a[10, 7] = convert(T, 15571226634087127616 // 774466927638876610083)
    a[10, 8] = convert(T, -133736375367792139885 // 4717207650164066625051)
    a[10, 9] = convert(T, 7461389216 // 501451974639)
    a[11, 1] = convert(T, 501256914705531962342417557181 // 14270506505142656332600844507392)
    a[11, 3] = convert(T, -1143766215625 // 132752960853408)
    a[11, 4] = convert(T, -6864570325 // 1185294293334)
    a[11, 5] = convert(T, 194348369382310456605879163404183 // 99893545535998594328205911551744)
    a[11, 6] = convert(T, -94634958447010580589908066176109375 // 27549212808177898050085930321520256)
    a[11, 7] = convert(T, -17006472665356285286219618514 // 155584463413110817059022733377)
    a[11, 8] = convert(T, 33530528814694461893884349656345 // 14270506505142656332600844507392)
    a[11, 9] = convert(T, -13439782155791134368 // 17777268379678341919)
    a[11, 10] = convert(T, 1441341768767571 // 13159456712985856)
    a[12, 1] = convert(
        T,
        parse(BigInt, "105854110734231079069010159870911189747853") //
            parse(BigInt, "5156624149476760916008179453333467046288864")
    )
    a[12, 3] = convert(T, -144579793509250000 // 19842290513127000261)
    a[12, 4] = convert(T, -101935644099967250 // 48188419817594143491)
    a[12, 5] = convert(
        T,
        parse(BigInt, "1585474394319811696785932424388196965") //
            parse(BigInt, "1709257457318830856936350991091849456")
    )
    a[12, 6] = convert(
        T,
        parse(BigInt, "-843499776333774172853009613469456309715703125") //
            parse(BigInt, "510505790798199330684809765880013237582597536")
    )
    a[12, 7] = convert(
        T,
        parse(BigInt, "-15057703799298260121553794369056896088480") //
            parse(BigInt, "714327132646734138085088291809720015274157")
    )
    a[12, 8] = convert(
        T,
        parse(BigInt, "1749840442221344572962864758990584360232600") //
            parse(BigInt, "1450300542040339007627300471250037606768743")
    )
    a[12, 9] = convert(T, -11255775246405733991656178432768 // 27206626483067760480757659602193)
    a[12, 10] = convert(T, 669010348769579696 // 7368057640845834597)
    a[12, 11] = convert(T, 4598083098752 // 858563707934367)
    a[13, 1] = convert(
        T,
        parse(BigInt, "-1639758773684715326849438048667467886824967397") //
            parse(BigInt, "11447568726280607813664651120965112496134881280")
    )
    a[13, 3] = convert(T, 3942453384375 // 314673684985856)
    a[13, 4] = convert(T, 11737114158175 // 1719466921529856)
    a[13, 5] = convert(
        T,
        -23710715033675876683332701739887457 //
            4940189888325748664958546898558976
    )
    a[13, 6] = convert(
        T,
        parse(BigInt, "498150575499633273684774666731162498301909124515625") //
            parse(BigInt, "87415924307623977386706008889913792042985180430336")
    )
    a[13, 7] = convert(
        T,
        parse(BigInt, "64881557768202140428371179540010005713998551") //
            parse(BigInt, "85896810580242200654071863296887242202224768")
    )
    a[13, 8] = convert(
        T,
        parse(BigInt, "-2336309182318568698279006266321563486172654055") //
            parse(BigInt, "18316109962048972501863441793544179993815810048")
    )
    a[13, 9] = convert(
        T,
        -493399374030747471036018890494175 //
            251658285736841065236836942273664
    )
    a[13, 10] = convert(T, 418285003077108927126515545155 // 455369916679568501838710898688)
    a[13, 11] = convert(T, -15171723902781457 // 63532954684873728)
    a[13, 12] = convert(T, 1501203688494867 // 9434957026426880)
    a[14, 1] = convert(
        T,
        parse(BigInt, "34188549803371802849576690267872548602326398788953") //
            parse(BigInt, "42496542183406636759747616530102745233754251202880")
    )
    a[14, 3] = convert(T, -18971246281693750 // 1138830954584356089)
    a[14, 4] = convert(T, -59230464334542700 // 2765732318276293359)
    a[14, 5] = convert(
        T,
        parse(BigInt, "5147939981309774383134903239728881770043") //
            parse(BigInt, "305929030949718561059100251282184099064")
    )
    a[14, 6] = convert(
        T,
        parse(BigInt, "-3625720213550267723370658302114678215563058405229078120") //
            parse(BigInt, "324512095420929759624784749347170583153994213035432256")
    )
    a[14, 7] = convert(
        T,
        parse(BigInt, "-60305503318319653518547439098565661266182518307816") //
            parse(BigInt, "17856872599361492097414471889911176856851308259643")
    )
    a[14, 8] = convert(
        T,
        parse(BigInt, "-1036461878759982363277481306266144563833492657780645") //
            parse(BigInt, "67994467493450618815596186448164392374006801924608")
    )
    a[14, 9] = convert(
        T,
        parse(BigInt, "128398681100219349205889126776607047000") //
            parse(BigInt, "7473801441221286756994805323613917077")
    )
    a[14, 10] = convert(T, -49156374556350058671822606102117 // 9039888303968618912866414995904)
    a[14, 11] = convert(T, 12253036339964386945 // 8828680926314891943)
    a[14, 12] = convert(T, -647188390508758231059 // 1092148506009694282240)
    a[14, 13] = convert(T, 10915833599872 // 368729913707897)
    a[15, 1] = convert(
        T,
        parse(BigInt, "-4939337286263213195547765488387521892799075623007291241961609516532") //
            parse(BigInt, "5408250052307451520718178852915698257207815452080611897685945761264")
    )
    a[15, 3] = convert(
        T,
        7588799849596321243074032368290625 //
            parse(BigInt, "3147217749590114939838670370597819616")
    )
    a[15, 4] = convert(
        T,
        16870665568420512953501332587233725 //
            955405388268427749593882076788623812
    )
    a[15, 5] = convert(
        T,
        parse(BigInt, "-808642515918378014850308582271476014669568437579087796060") //
            parse(BigInt, "54447992506702009927986632715967769032585338753056786562")
    )
    a[15, 6] = convert(
        T,
        parse(BigInt, "4610328329649866588704236006423149172472141907645890762410296050212") //
            parse(BigInt, "2135428689710103309390449198881479603148467934048051598947383737508")
    )
    a[15, 7] = convert(
        T,
        parse(BigInt, "4159963831215576225909381034291748993887819834160487158570788681") //
            parse(BigInt, "1040533184037697645660563795162185415624171583014576682740416336")
    )
    a[15, 8] = convert(
        T,
        parse(BigInt, "7381392142124351279433801934148706553542137071890521365664606664449580") //
            parse(BigInt, "259596002510757672994472584939953516345975141699869371088925396540699")
    )
    a[15, 9] = convert(
        T,
        parse(BigInt, "-3336834334584052813468828675971359774694437229547862706920") //
            parse(BigInt, "132102862435303266640535426836147775872819092781208127980")
    )
    a[15, 10] = convert(
        T,
        parse(BigInt, "426619379967412086875039012957475466130081426048213491790") //
            parse(BigInt, "55162410119399855550108207148248549410926885937244965785")
    )
    a[15, 11] = convert(
        T,
        parse(BigInt, "-630755628691078947314733435975762542732598947") //
            parse(BigInt, "333503232300511886435069380727586592765317456")
    )
    a[15, 12] = convert(
        T,
        parse(BigInt, "1522350657470125698997653827133798314909646891") //
            parse(BigInt, "1520094067152619944607524353149267399623188480")
    )
    a[15, 13] = convert(
        T,
        305575414262755427083262606101825880 //
            parse(BigInt, "65839748482572312891297405431209259829")
    )
    a[15, 14] = convert(
        T,
        parse(BigInt, "256624643108055110568255672032710477795") //
            parse(BigInt, "22874609758516552135947898572671559986304")
    )
    a[16, 1] = convert(
        T,
        parse(BigInt, "-571597862947184314270186718640978947715678864684269066846") //
            parse(BigInt, "2077055064880303907616135969012720011907767004397744786340")
    )
    a[16, 3] = convert(T, 66981514290625 // 1829501741761029)
    a[16, 4] = convert(T, 43495576635800 // 4443075658562499)
    a[16, 5] = convert(
        T,
        -127865248353371207265315478623656127 //
            10401415428935853634424440540325344
    )
    a[16, 6] = convert(
        T,
        parse(BigInt, "1316565142658075739557231574080234814338066993483960326560") //
            parse(BigInt, "92668695535091962564795912774190176478892159517481612467")
    )
    a[16, 7] = convert(
        T,
        parse(BigInt, "3881494143728609118531066904799685950051960514138645179820") //
            parse(BigInt, "2446349095978358868919950548516272963929118212742344026549")
    )
    a[16, 8] = convert(
        T,
        parse(BigInt, "162922667049680755852592453758428194006198229544701786842910") //
            parse(BigInt, "66288722243155885736983218667976563740242178853010092663614")
    )
    a[16, 9] = convert(
        T,
        parse(BigInt, "-43986024977384568043684084266385512680544563954") //
            parse(BigInt, "4922783599524658241955780540171948284522386185")
    )
    a[16, 10] = convert(
        T,
        parse(BigInt, "285912200202585226675651763671663063668290787") //
            parse(BigInt, "65371192072964016939690070594254881767827200")
    )
    a[16, 11] = convert(T, -6776815256667778089672518929 // 3693654613173093729492918708)
    a[16, 12] = convert(T, 398946554885847045598775476868169 // 344154261237450078839899047372800)
    a[16, 13] = convert(T, -76630698033396272 // 4432017119727044925)
    a[16, 14] = convert(T, 28401702316003037 // 1469612686944417840)
    a[16, 15] = convert(
        T,
        66049942462586341419969330578128801 //
            parse(BigInt, "12691068622536592094919763114637498325")
    )
    a[17, 1] = convert(
        T,
        parse(BigInt, "83940754497395557520874219603241359529066454343054832302344735") //
            parse(BigInt, "64192596456995578553872477759926464976144474354415663868673233")
    )
    a[17, 3] = convert(T, 892543892035485503125 // 51401651664490002607536)
    a[17, 4] = convert(T, -12732238157949399705325 // 686579204375687891972088)
    a[17, 5] = convert(
        T,
        parse(BigInt, "5290376174838819557032232941734928484252549") //
            parse(BigInt, "357179779572898187570048915214361602000384")
    )
    a[17, 6] = convert(
        T,
        parse(BigInt, "26873229338017506937199991804717456666650215387938173031932210") //
            parse(BigInt, "2863980005760296740624015421425947092438943496681472214589916")
    )
    a[17, 7] = convert(
        T,
        parse(BigInt, "-1976497866818803305857417297961598735637414137241493515492778650") //
            parse(BigInt, "378029217824623393200881653405474359138017953416246216408422692")
    )
    a[17, 8] = convert(
        T,
        parse(BigInt, "-1002860756304839757040188283199900676042073362417943601440986856950") //
            parse(BigInt, "20486915674765670626893195919603679319429068544972409068469849579")
    )
    a[17, 9] = convert(
        T,
        parse(BigInt, "87398661196965758104117684348440686081062878816711392590") //
            parse(BigInt, "2282122412587168891929052689609009868137678763277087160")
    )
    a[17, 10] = convert(
        T,
        parse(BigInt, "-7922242431969626895355493632206885458496418610471389") //
            parse(BigInt, "748272134517487495468365669337985635214015258726400")
    )
    a[17, 11] = convert(
        T,
        parse(BigInt, "2777643183645212014464950387658055285") //
            parse(BigInt, "1141545470045611737197667093465955392")
    )
    a[17, 12] = convert(
        T,
        parse(BigInt, "-1372659703515496442825084239977218110461") //
            parse(BigInt, "1313121960368535725613950174847107891200")
    )
    a[17, 13] = convert(T, 6144417902699179309851023 // 85608793932459282773805825)
    a[17, 14] = convert(T, 140294243355138853053241 // 64884622846351585391642880)
    a[17, 15] = convert(
        T,
        parse(BigInt, "168671028523891369934964082754523881107337") //
            parse(BigInt, "24062875279623260368388427013982199424119600")
    )
    b = [
        convert(T, 63818747 // 5262156900), zero(T), zero(T), zero(T), zero(T), zero(T),
        convert(T, 22555300000000 // 261366897038247),
        convert(T, 1696514453125 // 6717619827072),
        convert(T, -45359872 // 229764843),
        convert(T, 19174962087 // 94371046000),
        convert(T, -19310468 // 929468925),
        convert(T, 16089185487681 // 146694672924800),
        convert(T, 1592709632 // 41841694125),
        convert(T, 52675701958271 // 4527711056573100),
        convert(
            T,
            parse(BigInt, "12540904472870916741199505796420811396") //
                parse(BigInt, "2692319557780977037279406889319526430375")
        ),
        zero(T), zero(T),
    ]
    bp = [
        convert(T, 63818747 // 5262156900), zero(T), zero(T), zero(T), zero(T), zero(T),
        convert(T, 451106000000000 // 4965971043726693),
        convert(T, 8482572265625 // 26870479308288),
        convert(T, -181439488 // 689294529),
        convert(T, 57524886261 // 188742092000),
        convert(T, -38620936 // 929468925),
        convert(T, 144802669389129 // 586778691699200),
        convert(T, 6370838528 // 41841694125),
        convert(T, 368729913707897 // 4527711056573100),
        convert(
            T,
            parse(BigInt, "111940113324845802831946788738852162520696") //
                parse(BigInt, "1316544263754897771229629968877248424453375")
        ),
        convert(T, -113178587 // 12362232960),
        convert(T, 1 // 40),
    ]
    btilde = [
        convert(T, 63818747 // 5262156900 - 27121957 // 1594593000),
        zero(T), zero(T), zero(T), zero(T), zero(T),
        convert(T, 22555300000000 // 261366897038247 - 4006163300000 // 55441463008113),
        convert(T, 1696514453125 // 6717619827072 - 9466403125 // 25445529648),
        convert(T, -45359872 // 229764843 + 163199648 // 406149975),
        convert(T, 19174962087 // 94371046000 - 23359833 // 69636250),
        convert(T, -19310468 // 929468925 + 18491714 // 140828625),
        convert(T, 16089185487681 // 146694672924800 - 11052304606701 // 58344472186000),
        convert(T, 1592709632 // 41841694125 - 1191129152 // 44377554375),
        convert(T, 52675701958271 // 4527711056573100 - 2033811086741 // 124730332137000),
        convert(
            T,
            parse(BigInt, "12540904472870916741199505796420811396") //
                parse(BigInt, "2692319557780977037279406889319526430375") -
                parse(BigInt, "3616943474975740389660406409450169802") //
                parse(BigInt, "951830146690244407118982233597812374375")
        ),
        zero(T), zero(T),
    ]
    bptilde = [
        convert(T, 63818747 // 5262156900 - 27121957 // 1594593000),
        zero(T), zero(T), zero(T), zero(T), zero(T),
        convert(T, 451106000000000 // 4965971043726693 - 4217014000000 // 55441463008113),
        convert(T, 8482572265625 // 26870479308288 - 47332015625 // 101782118592),
        convert(T, -181439488 // 689294529 + 652798592 // 1218449925),
        convert(T, 57524886261 // 188742092000 - 70079499 // 139272500),
        convert(T, -38620936 // 929468925 + 36983428 // 140828625),
        convert(T, 144802669389129 // 586778691699200 - 99470741460309 // 233377888744000),
        convert(T, 6370838528 // 41841694125 - 4764516608 // 44377554375),
        convert(T, 368729913707897 // 4527711056573100 - 14236677607187 // 124730332137000),
        convert(
            T,
            parse(BigInt, "111940113324845802831946788738852162520696") //
                parse(BigInt, "1316544263754897771229629968877248424453375") -
                parse(BigInt, "198066487470143918516004831967805004004") //
                parse(BigInt, "2855490440070733221356946700793437123125")
        ),
        convert(T, -113178587 // 12362232960 - 1 // 50),
        convert(T, 1 // 40),
    ]
    c = T2[
        convert(T2, 1 // 50), convert(T2, 1 // 25), convert(T2, 1 // 10),
        convert(T2, 2 // 15), convert(T2, 4 // 25), convert(T2, 1 // 20),
        convert(T2, 1 // 5), convert(T2, 1 // 4), convert(T2, 1 // 3),
        convert(T2, 1 // 2), convert(T2, 5 // 9), convert(T2, 3 // 4),
        convert(T2, 6 // 7), convert(T2, 8437 // 8926), one(T2), one(T2),
    ]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function ERKN4Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 4, 4)
    a[2, 1] = convert(T, 1 // 32)
    a[3, 1] = convert(T, 19 // 600); a[3, 2] = convert(T, 16 // 75)
    a[4, 1] = convert(T, 32 // 315); a[4, 2] = convert(T, 58 // 315); a[4, 3] = convert(T, 3 // 14)
    b = [convert(T, 1 // 21), convert(T, 28 // 81), convert(T, 50 // 567), convert(T, 1 // 54)]
    bp = [convert(T, 1 // 14), convert(T, 32 // 81), convert(T, 250 // 567), convert(T, 5 // 54)]
    btilde = [
        convert(T, 1 // 21 - 14 // 375), convert(T, 28 // 81 - 136 // 375),
        convert(T, 50 // 567 - 2 // 25), convert(T, 1 // 54 - 1 // 50),
    ]
    bptilde = [
        convert(T, 1 // 14 - 17 // 231), convert(T, 32 // 81 - 116 // 297),
        convert(T, 250 // 567 - 925 // 2079), convert(T, 5 // 54 - 1 // 11),
    ]
    c = T2[convert(T2, 1 // 4), convert(T2, 7 // 10), one(T2)]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function ERKN5Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 4, 4)
    a[2, 1] = convert(T, 1 // 8)
    a[3, 1] = convert(T, 2907 // 343000); a[3, 2] = convert(T, 1216 // 42875)
    a[4, 1] = convert(T, 6624772 // Int64(128538819))
    a[4, 2] = convert(T, 6273905 // Int64(54121608))
    a[4, 3] = convert(T, Int64(210498365) // Int64(1028310552))
    b = [
        convert(T, 479 // 5016), convert(T, 235 // 1776),
        convert(T, 145775 // 641744), convert(T, 309519 // 6873416),
    ]
    bp = [
        convert(T, 479 // 5016), convert(T, 235 // 888),
        convert(T, 300125 // 962616), convert(T, 2255067 // 6873416),
    ]
    btilde = [
        convert(T, 479 // 5016 - 184883 // 2021250),
        convert(T, 235 // 1776 - 411163 // 3399375),
        convert(T, 145775 // 641744 - 6 // 25),
        convert(T, 309519 // 6873416 - 593028 // Int64(12464375)),
    ]
    bptilde = T[]
    c = T2[convert(T2, 1 // 2), convert(T2, 19 // 70), convert(T2, 44 // 51)]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, true, 2)
end

function ERKN7Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 7, 7)
    a[2, 1] = convert(T, 5107771 // 767472028)
    a[3, 1] = convert(T, 5107771 // 575604021); a[3, 2] = convert(T, 16661485 // 938806552)
    a[4, 1] = convert(T, 325996677 // 876867260); a[4, 2] = convert(T, -397622579 // 499461366)
    a[4, 3] = convert(T, 541212017 // 762248206)
    a[5, 1] = convert(T, 82243160 // 364375691); a[5, 2] = convert(T, -515873404 // 1213273815)
    a[5, 3] = convert(T, 820109726 // 1294837243); a[5, 4] = convert(T, 36245507 // 242779260)
    a[6, 1] = convert(T, 3579594 // 351273191); a[6, 2] = convert(T, 34292133 // 461028419)
    a[6, 3] = convert(T, 267156948 // 2671391749); a[6, 4] = convert(T, 22665163 // 1338599875)
    a[6, 5] = convert(T, -3836509 // 1614789462)
    a[7, 1] = convert(T, 53103334 // 780726093); a[7, 3] = convert(T, 352190060 // 1283966121)
    a[7, 4] = convert(T, 37088117 // 2206150964); a[7, 5] = convert(T, 7183323 // 1828127386)
    a[7, 6] = convert(T, 187705681 // 1370684829)
    b = [
        convert(T, 53103334 // 780726093), zero(T), convert(T, 352190060 // 1283966121),
        convert(T, 37088117 // 2206150964), convert(T, 7183323 // 1828127386),
        convert(T, 187705681 // 1370684829), zero(T),
    ]
    bp = [
        convert(T, 53103334 // 780726093), zero(T), convert(T, 244481296 // 685635505),
        convert(T, 41493456 // 602487871), convert(T, -45498718 // 926142189),
        convert(T, 1625563237 // 4379140271), convert(T, 191595797 // 1038702495),
    ]
    btilde = [
        convert(T, 53103334 // 780726093 - 41808761 // 935030896), zero(T),
        convert(T, 352190060 // 1283966121 - 46261019 // 135447428),
        convert(T, 37088117 // 2206150964 - 289298425 // 1527932372),
        convert(T, 7183323 // 1828127386 + 52260067 // 3104571287),
        convert(T, 187705681 // 1370684829 + 49872919 // 848719175), zero(T),
    ]
    bptilde = [
        convert(T, 53103334 // 780726093 - 41808761 // 935030896), zero(T),
        convert(T, 244481296 // 685635505 - 224724272 // 506147085),
        convert(T, 41493456 // 602487871 - 2995752066 // 3862177123),
        convert(T, -45498718 // 926142189 - 170795979 // 811534085),
        convert(T, 1625563237 // 4379140271 + 177906423 // 1116903503),
        convert(T, 191595797 // 1038702495 + 655510901 // 2077404990),
    ]
    c = T2[
        convert(T2, 108816483 // 943181462), convert(T2, 108816483 // 471590731),
        convert(T2, 151401202 // 200292705), convert(T2, 682035803 // 631524599),
        convert(T2, 493263404 // 781610081), one(T2),
    ]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function Nystrom5VelocityIndependentTableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 4, 4)
    a[2, 1] = convert(T, 1 // 50)
    a[3, 1] = convert(T, -1 // 27); a[3, 2] = convert(T, 7 // 27)
    a[4, 1] = convert(T, 3 // 10); a[4, 2] = convert(T, -2 // 35); a[4, 3] = convert(T, 9 // 35)
    b = [convert(T, 14 // 336), convert(T, 100 // 336), convert(T, 54 // 336), zero(T)]
    bp = [convert(T, 14 // 336), convert(T, 125 // 336), convert(T, 162 // 336), convert(T, 35 // 336)]
    btilde = T[]
    bptilde = T[]
    c = T2[convert(T2, 1 // 5), convert(T2, 2 // 3), one(T2)]
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

"""
    NystromVDTableau{T, T2}

Tableau for velocity-dependent Nyström methods.
Fields:
- `a`: nstages × nstages lower-triangular position coupling matrix
- `abar`: nstages × nstages lower-triangular velocity coupling matrix
- `b`: position update weights (length nstages)
- `bp`: velocity update weights (length nstages)
- `btilde`: embedded position error weights (empty if non-adaptive)
- `bptilde`: embedded velocity error weights (empty if non-adaptive)
- `c`: time nodes for stages 2..nstages (length nstages-1); c[i] is node for stage i+1
- `nf_per_step`: number of f1 evaluations to count per step (default = nstages, i.e.,
  loop stages k2..kN plus fsallast; set to nstages-1 to exclude fsallast from the count)
"""
struct NystromVDTableau{T, T2}
    a::Matrix{T}
    abar::Matrix{T}
    b::Vector{T}
    bp::Vector{T}
    btilde::Vector{T}
    bptilde::Vector{T}
    c::Vector{T2}
    nf_per_step::Int
end

function Nystrom4VelocityIndependentTableau(::Type{T}, ::Type{T2}) where {T, T2}
    # 3 stages, velocity-independent: kᵢ = f1(duprev, kuᵢ, p, t+cᵢ*dt)
    # Coefficients from perform_step!:
    #   c = [1/2, 1]; a[2,1]=1/8, a[3,2]=1/2
    #   b = [1/6, 2/6, 0], bp = [1/6, 4/6, 1/6]
    nstages = 3
    a = zeros(T, nstages, nstages)
    a[2, 1] = convert(T, 1 // 8)
    a[3, 2] = convert(T, 1 // 2)
    b = [convert(T, 1 // 6), convert(T, 2 // 6), zero(T)]
    bp = [convert(T, 1 // 6), convert(T, 4 // 6), convert(T, 1 // 6)]
    btilde = T[]  # non-adaptive
    bptilde = T[]
    c = [convert(T2, 1 // 2), one(T2)]  # c for stages 2,3
    return NystromVITableau(a, b, bp, btilde, bptilde, c, false, 2)
end

function RKN4Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    # 3 stages, velocity-dependent
    # k2 at c[1]=1/2: ku = uprev + dt*(1/2)*duprev + dt²*(1/8)*k1
    #                 kdu = duprev + dt*(1/2)*k1
    # k3 at c[2]=1:   ku = uprev + dt*1*duprev + dt²*(1/2)*k2
    #                 kdu = duprev + dt*1*k2
    # u = uprev + dt*duprev + dt²*(1/6*k1 + 2/6*k2 + 0*k3)
    # du = duprev + dt*(1/6*k1 + 4/6*k2 + 1/6*k3)
    # nf_per_step=2: matches original counting convention (k2+k3 only, not fsallast)
    nstages = 3
    a = zeros(T, nstages, nstages)
    a[2, 1] = convert(T, 1 // 8)
    a[3, 2] = convert(T, 1 // 2)
    abar = zeros(T, nstages, nstages)
    abar[2, 1] = convert(T, 1 // 2)
    abar[3, 2] = convert(T, 1 // 1)
    b = [convert(T, 1 // 6), convert(T, 2 // 6), zero(T)]
    bp = [convert(T, 1 // 6), convert(T, 4 // 6), convert(T, 1 // 6)]
    btilde = T[]  # non-adaptive
    bptilde = T[]
    c = [convert(T2, 1 // 2), one(T2)]  # c for stages 2,3
    return NystromVDTableau(a, abar, b, bp, btilde, bptilde, c, nstages - 1)
end

function Nystrom4Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    # 4 stages, velocity-dependent
    # From perform_step! Nystrom4ConstantCache:
    # k2 at c[1]=1/2: ku = uprev + dt*(1/2)*duprev + dt²*(1/8)*k1
    #                 kdu = duprev + dt*(1/2)*k1
    # k3 at c[2]=1/2: ku = uprev + dt*(1/2)*duprev + dt²*(1/8)*k1  (same ku as k2!)
    #                 kdu = duprev + dt*(1/2)*k2
    # k4 at c[3]=1:   ku = uprev + dt*1*duprev + dt²*(1/2)*k3
    #                 kdu = duprev + dt*1*k3
    # u = uprev + dt*duprev + dt²*(1/6*(k1+k2+k3))   [b4=0]
    # du = duprev + dt*(1/6*(k1+k4) + 2/6*(k2+k3))
    nstages = 4
    a = zeros(T, nstages, nstages)
    a[2, 1] = convert(T, 1 // 8)
    a[3, 1] = convert(T, 1 // 8)   # a[3,2] = 0
    a[4, 3] = convert(T, 1 // 2)
    abar = zeros(T, nstages, nstages)
    abar[2, 1] = convert(T, 1 // 2)
    abar[3, 2] = convert(T, 1 // 2)
    abar[4, 3] = convert(T, 1 // 1)
    b = [convert(T, 1 // 6), convert(T, 1 // 6), convert(T, 1 // 6), zero(T)]
    bp = [convert(T, 1 // 6), convert(T, 2 // 6), convert(T, 2 // 6), convert(T, 1 // 6)]
    btilde = T[]  # non-adaptive
    bptilde = T[]
    c = [convert(T2, 1 // 2), convert(T2, 1 // 2), one(T2)]  # c for stages 2,3,4
    return NystromVDTableau(a, abar, b, bp, btilde, bptilde, c, nstages)
end

function FineRKN4Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 5, 5)
    a[2, 1] = convert(T, 2 // 81)
    a[3, 1] = convert(T, 1 // 36); a[3, 2] = convert(T, 1 // 36)
    a[4, 1] = convert(T, 9 // 128); a[4, 3] = convert(T, 27 // 128)
    a[5, 1] = convert(T, 11 // 60); a[5, 2] = convert(T, -3 // 20)
    a[5, 3] = convert(T, 9 // 25); a[5, 4] = convert(T, 8 // 75)
    abar = zeros(T, 5, 5)
    abar[2, 1] = convert(T, 2 // 9)
    abar[3, 1] = convert(T, 1 // 12); abar[3, 2] = convert(T, 1 // 4)
    abar[4, 1] = convert(T, 69 // 128); abar[4, 2] = convert(T, -243 // 128); abar[4, 3] = convert(T, 135 // 64)
    abar[5, 1] = convert(T, -17 // 12); abar[5, 2] = convert(T, 27 // 4)
    abar[5, 3] = convert(T, -27 // 5); abar[5, 4] = convert(T, 16 // 15)
    b = [convert(T, 19 // 180), zero(T), convert(T, 63 // 200), convert(T, 16 // 225), convert(T, 1 // 120)]
    bp = [convert(T, 1 // 9), zero(T), convert(T, 9 // 20), convert(T, 16 // 45), convert(T, 1 // 12)]
    btilde = [
        convert(T, 25 // 1116), zero(T), convert(T, -63 // 1240),
        convert(T, 64 // 1395), convert(T, -13 // 744),
    ]
    bptilde = [
        convert(T, 2 // 125), zero(T), convert(T, -27 // 625),
        convert(T, 32 // 625), convert(T, -3 // 125),
    ]
    c = T2[convert(T2, 2 // 9), convert(T2, 1 // 3), convert(T2, 3 // 4), one(T2)]
    return NystromVDTableau(a, abar, b, bp, btilde, bptilde, c, 5)
end

function FineRKN5Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    a = zeros(T, 7, 7)
    a[2, 1] = convert(T, 32 // 1521)
    a[3, 1] = convert(T, 4 // 169); a[3, 2] = convert(T, 4 // 169)
    a[4, 1] = convert(T, 175 // 5184); a[4, 3] = convert(T, 1625 // 5184)
    a[5, 1] = convert(T, -342497279 // 5618900760); a[5, 2] = convert(T, 6827067 // 46824173)
    a[5, 3] = convert(T, 35048741 // 102161832); a[5, 4] = convert(T, -2201514 // 234120865)
    a[6, 1] = convert(T, -7079 // 52152); a[6, 2] = convert(T, 767 // 2173)
    a[6, 3] = convert(T, 14027 // 52152); a[6, 4] = convert(T, 30 // 2173)
    a[7, 1] = convert(T, 4817 // 51600); a[7, 3] = convert(T, 388869 // 1216880)
    a[7, 4] = convert(T, 3276 // 23575); a[7, 5] = convert(T, -1142053 // 22015140)
    abar = zeros(T, 7, 7)
    abar[2, 1] = convert(T, 8 // 39)
    abar[3, 1] = convert(T, 1 // 13); abar[3, 2] = convert(T, 3 // 13)
    abar[4, 1] = convert(T, 7385 // 6912); abar[4, 2] = convert(T, -9425 // 2304)
    abar[4, 3] = convert(T, 13325 // 3456)
    abar[5, 1] = convert(T, 223324757 // 91364240); abar[5, 2] = convert(T, -174255393 // 18272848)
    abar[5, 3] = convert(T, 382840094 // 46824173); abar[5, 4] = convert(T, -39627252 // 234120865)
    abar[6, 1] = convert(T, 108475 // 36464); abar[6, 2] = convert(T, -9633 // 848)
    abar[6, 3] = convert(T, 7624604 // 806183); abar[6, 4] = convert(T, 8100 // 49979)
    abar[6, 5] = convert(T, -4568212 // 19446707)
    abar[7, 1] = convert(T, 4817 // 51600); abar[7, 3] = convert(T, 1685099 // 3650640)
    abar[7, 4] = convert(T, 19656 // 23575); abar[7, 5] = convert(T, -53676491 // 88060560)
    abar[7, 6] = convert(T, 53 // 240)
    b = [
        convert(T, 4817 // 51600), zero(T), convert(T, 388869 // 1216880),
        convert(T, 3276 // 23575), convert(T, -1142053 // 22015140), zero(T), zero(T),
    ]
    bp = [
        convert(T, 4817 // 51600), zero(T), convert(T, 1685099 // 3650640),
        convert(T, 19656 // 23575), convert(T, -53676491 // 88060560),
        convert(T, 53 // 240), zero(T),
    ]
    btilde = [
        convert(T, 8151 // 2633750), zero(T), convert(T, -1377519 // 186334750),
        convert(T, 586872 // 28879375), convert(T, -36011118 // 2247378875), zero(T), zero(T),
    ]
    bptilde = [
        convert(T, 8151 // 2633750), zero(T), convert(T, -5969249 // 559004250),
        convert(T, 3521232 // 28879375), convert(T, -846261273 // 4494757750),
        convert(T, 4187 // 36750), convert(T, -1 // 25),
    ]
    c = T2[
        convert(T2, 8 // 39), convert(T2, 4 // 13), convert(T2, 5 // 6),
        convert(T2, 43 // 47), one(T2), one(T2),
    ]
    return NystromVDTableau(a, abar, b, bp, btilde, bptilde, c, 7)
end

struct IRKN3ConstantCache{T, T2} <: NystromConstantCache
    bconst1::T
    bconst2::T
    c1::T2
    a21::T
    b1::T
    b2::T
    bbar1::T
    bbar2::T
end

function IRKN3ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    bconst1 = convert(T, 1.5)
    bconst2 = convert(T, -0.5)
    c1 = convert(T2, 0.5)
    a21 = convert(T, 0.125)
    b1 = convert(T, 0.6666666666666666)
    b2 = convert(T, 0.8333333333333334)
    bbar1 = convert(T, 0.3333333333333333)
    bbar2 = convert(T, 0.4166666666666667)
    return IRKN3ConstantCache(bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2)
end

function IRKN3ConstantCache(T::Type, T2::Type)
    bconst1 = convert(T, 3 // 2)
    bconst2 = convert(T, -1 // 2)
    c1 = convert(T2, 1 // 2)
    a21 = convert(T, 1 // 8)
    b1 = convert(T, 2 // 3)
    b2 = convert(T, 5 // 6)
    bbar1 = convert(T, 1 // 3)
    bbar2 = convert(T, 5 // 12)
    return IRKN3ConstantCache(bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2)
end

struct IRKN4ConstantCache{T, T2} <: NystromConstantCache
    bconst1::T
    bconst2::T
    c1::T2
    c2::T2
    a21::T
    # a31::T
    a32::T
    b1::T
    b2::T
    b3::T
    bbar1::T
    bbar2::T
    bbar3::T
end

function IRKN4ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    bconst1 = convert(T, 1.5)
    bconst2 = convert(T, -0.5)
    c1 = convert(T2, 0.25)
    c2 = convert(T2, 0.75)
    a21 = convert(T, 0.03125)
    # a31     = convert(T,0)
    a32 = convert(T, 0.28125)
    b1 = convert(T, 1.0555555555555556)
    b2 = convert(T, -0.16666666666666666)
    b3 = convert(T, 0.6111111111111112)
    bbar1 = convert(T, -0.05555555555555555)
    bbar2 = convert(T, 0.2916666666666667)
    bbar3 = convert(T, 0.125)
    return IRKN4ConstantCache(bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3)
end

function IRKN4ConstantCache(T::Type, T2::Type)
    bconst1 = convert(T, 3 // 2)
    bconst2 = convert(T, -1 // 2)
    c1 = convert(T2, 1 // 4)
    c2 = convert(T2, 3 // 4)
    a21 = convert(T, 1 // 32)
    # a31     = convert(T,0)
    a32 = convert(T, 9 // 32)
    b1 = convert(T, 19 // 18)
    b2 = convert(T, -1 // 6)
    b3 = convert(T, 11 // 18)
    bbar1 = convert(T, -1 // 18)
    bbar2 = convert(T, 7 // 24)
    bbar3 = convert(T, 1 // 8)
    return IRKN4ConstantCache(bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3)
end


# DPRKN6 carries the standard Nyström Butcher arrays (so it runs through the generic
# velocity-independent perform_step) plus its specialized "free" 6th-order dense-output
# coefficients r*/rp*, which the dedicated interpolant in interpolants.jl consumes.
struct DPRKN6Tableau{T, T2}
    a::Matrix{T}
    b::Vector{T}
    bp::Vector{T}
    btilde::Vector{T}
    bptilde::Vector{T}
    c::Vector{T2}
    pos_only_error::Bool
    # Dense-output polynomial coefficients. R[i, k+1] (resp. Rp) is the Θ^k coefficient
    # of bᵢ(Θ) (resp. bpᵢ(Θ)). Row 2 is unused (k₂ does not enter the interpolant).
    R::Matrix{T}
    Rp::Matrix{T}
    # See NystromVITableau.kshortsize. DPRKN6 uses 3 to cache stage derivatives for
    # its specialized 6th-order interpolant.
    kshortsize::Int
end

function DPRKN6Tableau(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.12929590313670442)
    c2 = convert(T2, 0.25859180627340883)
    c3 = convert(T2, 0.67029708261548)
    c4 = convert(T2, 0.9)
    c5 = convert(T2, 1.0)
    a21 = convert(T, 0.008358715283968025)
    a31 = convert(T, 0.011144953711957367)
    a32 = convert(T, 0.022289907423914734)
    a41 = convert(T, 0.1454747428010918)
    a42 = convert(T, -0.22986064052264749)
    a43 = convert(T, 0.3090349872029675)
    a51 = convert(T, -0.20766826295078997)
    a52 = convert(T, 0.6863667842925143)
    a53 = convert(T, -0.19954927787234925)
    a54 = convert(T, 0.12585075653062489)
    a61 = convert(T, 0.07811016144349478)
    a63 = convert(T, 0.2882917411897668)
    a64 = convert(T, 0.12242553717457041)
    a65 = convert(T, 0.011172560192168035)
    b1 = convert(T, 0.07811016144349478)
    b3 = convert(T, 0.2882917411897668)
    b4 = convert(T, 0.12242553717457041)
    b5 = convert(T, 0.011172560192168035)
    bp1 = convert(T, 0.07811016144349478)
    bp3 = convert(T, 0.3888434787059826)
    bp4 = convert(T, 0.3713207579288423)
    bp5 = convert(T, 0.11172560192168035)
    bp6 = convert(T, 0.05)
    btilde1 = convert(T, -0.9807490989269235)
    btilde2 = convert(T, 2.406751371924452)
    btilde3 = convert(T, -1.559600370364267)
    btilde4 = convert(T, 0.12242553717457041)
    btilde5 = convert(T, 0.011172560192168035)
    bptilde1 = convert(T, 0.023504273504273504)
    bptilde3 = convert(T, -0.07242330719764424)
    bptilde4 = convert(T, 0.17543989844952962)
    bptilde5 = convert(T, -0.2765208647561589)
    bptilde6 = convert(T, 0.15)
    r14 = convert(T, 0.21367521367521367)
    r13 = convert(T, -0.9066951566951567)
    r12 = convert(T, 1.5161443494776827)
    r11 = convert(T, -1.245014245014245)
    r10 = convert(T, 0.5)
    r34 = convert(T, -0.6583937017967658)
    r33 = convert(T, 2.5384011164109506)
    r32 = convert(T, -3.577652872294921)
    r31 = convert(T, 1.9859371988705032)
    r44 = convert(T, 1.5949081677229964)
    r43 = convert(T, -5.164133553908094)
    r42 = convert(T, 5.547586751052329)
    r41 = convert(T, -1.8559358276926614)
    r54 = convert(T, -2.513826043237808)
    r53 = convert(T, 7.273336685101391)
    r52 = convert(T, -6.926987319144182)
    r51 = convert(T, 2.178649237472767)
    r64 = convert(T, 1.3636363636363635)
    r63 = convert(T, -3.7409090909090907)
    r62 = convert(T, 3.440909090909091)
    r61 = convert(T, -1.0636363636363637)
    rp14 = convert(T, 1.2820512820512822)
    rp13 = convert(T, -4.533475783475783)
    rp12 = convert(T, 6.064577397910731)
    rp11 = convert(T, -3.735042735042735)
    rp10 = convert(T, 1)
    rp34 = convert(T, -3.950362210780595)
    rp33 = convert(T, 12.692005582054751)
    rp32 = convert(T, -14.310611489179683)
    rp31 = convert(T, 5.95781159661151)
    rp44 = convert(T, 9.56944900633798)
    rp43 = convert(T, -25.820667769540467)
    rp42 = convert(T, 22.190347004209315)
    rp41 = convert(T, -5.567807483077984)
    rp54 = convert(T, -15.082956259426847)
    rp53 = convert(T, 36.366683425506956)
    rp52 = convert(T, -27.707949276576727)
    rp51 = convert(T, 6.5359477124183005)
    rp64 = convert(T, 8.181818181818182)
    rp63 = convert(T, -18.704545454545453)
    rp62 = convert(T, 13.763636363636364)
    rp61 = convert(T, -3.190909090909091)
    return _assemble_dprkn6_tableau(
        T, T2, c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51,
        a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
        bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4,
        btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6,
        r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41,
        r54,
        r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10,
        rp34,
        rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51,
        rp64, rp63, rp62, rp61
    )
end

function DPRKN6Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    R = sqrt(big(8581))
    c1 = convert(T2, (209 - R) / 900)
    c2 = convert(T2, (209 - R) / 450)
    c3 = convert(T2, (209 + R) / 450)
    c4 = convert(T2, 9 // 10)
    c5 = convert(T2, 1)
    a21 = convert(T, (26131 - 209R) / 81_0000)
    a31 = convert(T, (26131 - 209R) / 60_7500)
    a32 = convert(T, (26131 - 209R) / 30_3750)
    a41 = convert(T, (980403512254 + 7781688431R) / 116944_6992_1875)
    a42 = convert(T, -(126288_4486208 + 153854_81287R) / 116944_6992_1875)
    a43 = convert(T, (7166_233_891_441 + 786_945_632_99R) / 46_777_879_687_500)
    a51 = convert(T, -9(329260 + 3181R) / 2704_0000)
    a52 = convert(T, 27(35129 + 3331R) / 1352_0000)
    a53 = convert(T, -27(554358343 + 31040327R) / 46406048_0000)
    a54 = convert(T, 153(8555_257 - 67973R) / 274592_0000)
    a61 = convert(T, 329 // 4212)
    # a62      = convert(T,0)
    a63 = convert(T, (8411_9543 + 366_727R) / 4096_22616)
    a64 = convert(T, (8411_9543 - 366_727R) / 4096_22616)
    a65 = convert(T, 200 // 17901)
    b1 = convert(T, 329 // 4212)
    # b2       = convert(T,0)
    b3 = a63
    b4 = a64
    b5 = convert(T, 200 // 17901)
    # b6       = convert(T,0)
    bp1 = b1
    # bp2      = b2
    bp3 = convert(T, (389225579 + 96856R) / 10_2405_6540)
    bp4 = convert(T, (389225579 - 96856R) / 10_2405_6540)
    bp5 = convert(T, 2000 // 17901)
    bp6 = convert(T, 1 // 20)
    btilde1 = convert(T, 329 // 4212 - (2701 + 23R) / 4563)
    btilde2 = convert(T, (9829 + 131R) / 9126)
    btilde3 = convert(T, (8411_9543 + 366_727R) / 4096_22616 - 5(1798 + 17R) / 9126)
    btilde4 = b4
    btilde5 = b5
    # btilde6  = convert(T,0)
    bptilde1 = convert(T, 329 // 4212 - 115 // 2106)
    # btildep2 = convert(T,0)
    bptilde3 = convert(
        T,
        (389225579 + 96856R) / 10_2405_6540 -
            (8411_9543 + 366_727R) / 2560_14135
    )
    bptilde4 = convert(
        T,
        (389225579 - 96856R) / 10_2405_6540 -
            (8411_9543 - 366_727R) / 2560_14135
    )
    bptilde5 = convert(T, 2000 // 17901 - 6950 // 17901)
    bptilde6 = convert(T, 1 // 20 + 1 // 10)
    r14 = convert(T, 900 // 4212)
    r13 = convert(T, -3819 // 4212)
    r12 = convert(T, 6386 // 4212)
    r11 = convert(T, -5244 // 4212)
    r10 = convert(T, 2106 // 4212)
    r34 = convert(T, 1800 * (5860823 - 152228R) / 22529243880)
    r33 = convert(T, -6 * (4929647204 - 156109769R) / 22529243880)
    r32 = convert(T, (22190560391 - 1109665151R) / 22529243880)
    r31 = convert(T, 18 * (81356461 + 25954829R) / 22529243880)
    r44 = convert(T, 1800 * (5860823 + 152228R) / 22529243880)
    r43 = convert(T, -6 * (4929647204 + 156109769R) / 22529243880)
    r42 = convert(T, (22190560391 + 1109665151R) / 22529243880)
    r41 = convert(T, 18 * (81356461 - 25954829R) / 22529243880)
    r54 = convert(T, -200 * 225 // 17901)
    r53 = convert(T, 200 * 651 // 17901)
    r52 = convert(T, -200 * 620 // 17901)
    r51 = convert(T, 200 * 195 // 17901)
    r64 = convert(T, 15 // 11)
    r63 = convert(T, -823 // 220)
    r62 = convert(T, 757 // 220)
    r61 = convert(T, -117 // 110)
    rp14 = convert(T, 5400 // 4212)
    rp13 = convert(T, -19095 // 4212)
    rp12 = convert(T, 25544 // 4212)
    rp11 = convert(T, -15732 // 4212)
    rp10 = convert(T, 1)
    rp34 = convert(T, 5400 * (5860823 - 152228R) / 11264621940)
    rp33 = convert(T, -15 * (4929647204 - 156109769R) / 11264621940)
    rp32 = convert(T, 2 * (22190560391 - 1109665151R) / 11264621940)
    rp31 = convert(T, 27 * (81356461 + 25954829R) / 11264621940)
    rp44 = convert(T, 5400 * (5860823 + 152228R) / 11264621940)
    rp43 = convert(T, -15 * (4929647204 + 156109769R) / 11264621940)
    rp42 = convert(T, 2 * (22190560391 + 1109665151R) / 11264621940)
    rp41 = convert(T, 27 * (81356461 - 25954829R) / 11264621940)
    rp54 = convert(T, -1000 * 270 // 17901)
    rp53 = convert(T, 1000 * 651 // 17901)
    rp52 = convert(T, -1000 * 496 // 17901)
    rp51 = convert(T, 1000 * 117 // 17901)
    rp64 = convert(T, 1800 // 220)
    rp63 = convert(T, -4115 // 220)
    rp62 = convert(T, 3028 // 220)
    rp61 = convert(T, -702 // 220)
    return _assemble_dprkn6_tableau(
        T, T2, c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51,
        a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
        bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4,
        btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6,
        r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41,
        r54,
        r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10,
        rp34,
        rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51,
        rp64, rp63, rp62, rp61
    )
end

# Assemble DPRKN6's scalar Butcher coefficients into the array form the generic
# velocity-independent perform_step consumes, packaged with the dense-output coeffs.
function _assemble_dprkn6_tableau(
        ::Type{T}, ::Type{T2}, c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51,
        a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
        bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4,
        btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6,
        r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41,
        r54, r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10,
        rp34, rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51,
        rp64, rp63, rp62, rp61
    ) where {T, T2}
    a = zeros(T, 6, 6)
    a[2, 1] = a21
    a[3, 1] = a31; a[3, 2] = a32
    a[4, 1] = a41; a[4, 2] = a42; a[4, 3] = a43
    a[5, 1] = a51; a[5, 2] = a52; a[5, 3] = a53; a[5, 4] = a54
    a[6, 1] = a61; a[6, 3] = a63; a[6, 4] = a64; a[6, 5] = a65
    b = T[b1, zero(T), b3, b4, b5, zero(T)]
    bp = T[bp1, zero(T), bp3, bp4, bp5, bp6]
    btilde = T[btilde1, btilde2, btilde3, btilde4, btilde5, zero(T)]
    bptilde = T[bptilde1, zero(T), bptilde3, bptilde4, bptilde5, bptilde6]
    c = T2[c1, c2, c3, c4, c5]
    # Pack dense-output coefficients into 6×5 matrices. Row 2 (k₂) is unused.
    # b₁(Θ) is a full degree-4 polynomial; the other bᵢ(Θ) have a zero constant term,
    # written explicitly here so the interpolant evaluates `@evalpoly` uniformly.
    z = zero(T)
    R = T[
        r10  r11  r12  r13  r14
        z    z    z    z    z
        z    r31  r32  r33  r34
        z    r41  r42  r43  r44
        z    r51  r52  r53  r54
        z    r61  r62  r63  r64
    ]
    Rp = T[
        rp10 rp11 rp12 rp13 rp14
        z    z    z    z    z
        z    rp31 rp32 rp33 rp34
        z    rp41 rp42 rp43 rp44
        z    rp51 rp52 rp53 rp54
        z    rp61 rp62 rp63 rp64
    ]
    return DPRKN6Tableau(a, b, bp, btilde, bptilde, c, false, R, Rp, 3)
end
