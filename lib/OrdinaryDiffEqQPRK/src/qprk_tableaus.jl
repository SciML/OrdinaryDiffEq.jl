struct QPRK98Tableau{T, T2}
    d2::T2
    d3::T2
    d4::T2
    d5::T2
    d6::T2
    d7::T2
    d8::T2
    d9::T2
    d10::T2
    d11::T2
    d12::T2
    d13::T2
    d14::T2
    ϵ1::T
    ϵ8::T
    ϵ9::T
    ϵ10::T
    ϵ11::T
    ϵ12::T
    ϵ13::T
    ϵ14::T
    ϵ15::T
    w1::T
    w8::T
    w9::T
    w10::T
    w11::T
    w12::T
    w13::T
    w14::T
    w15::T
    w16::T
    b21::T
    b31::T
    b32::T
    b41::T
    b43::T
    b51::T
    b53::T
    b54::T
    b61::T
    b64::T
    b65::T
    b71::T
    b74::T
    b75::T
    b76::T
    b81::T
    b86::T
    b87::T
    b91::T
    b96::T
    b97::T
    b98::T
    b10_1::T
    b10_6::T
    b10_7::T
    b10_8::T
    b10_9::T
    b11_1::T
    b11_6::T
    b11_7::T
    b11_8::T
    b11_9::T
    b11_10::T
    b12_1::T
    b12_6::T
    b12_7::T
    b12_8::T
    b12_9::T
    b12_10::T
    b12_11::T
    b13_1::T
    b13_6::T
    b13_7::T
    b13_8::T
    b13_9::T
    b13_10::T
    b13_11::T
    b13_12::T
    b14_1::T
    b14_6::T
    b14_7::T
    b14_8::T
    b14_9::T
    b14_10::T
    b14_11::T
    b14_12::T
    b14_13::T
    b15_1::T
    b15_6::T
    b15_7::T
    b15_8::T
    b15_9::T
    b15_10::T
    b15_11::T
    b15_12::T
    b15_13::T
    b15_14::T
    b16_1::T
    b16_6::T
    b16_7::T
    b16_8::T
    b16_9::T
    b16_10::T
    b16_11::T
    b16_12::T
    b16_13::T
    b16_14::T
end

@fold function QPRK98Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    d2 = convert(T, BigInt(3) // BigInt(184))
    d3 = convert(T, BigInt(1005991230711768) // BigInt(26022242251007185))
    d4 = convert(T, BigInt(1003144597305563) // BigInt(17299071752613603))
    d5 = convert(T, BigInt(3230) // BigInt(3269))
    d6 = convert(T, BigInt(3634874590315107) // BigInt(41460143603477948))
    d7 = convert(T, BigInt(840156096102871) // BigInt(4026814668841799))
    d8 = convert(T, BigInt(703) // BigInt(2847))
    d9 = convert(T, BigInt(703) // BigInt(3796))
    d10 = convert(T, BigInt(2135) // BigInt(9119))
    d11 = convert(T, BigInt(1531) // BigInt(2485))
    d12 = convert(T, BigInt(15434796193306477) // BigInt(18528760750836549))
    d13 = convert(T, BigInt(1999) // BigInt(2000))
    d14 = convert(T, BigInt(6594) // BigInt(6595))
    ϵ1 = convert(T2, BigInt(32020052252709) // BigInt(309017382031080920))
    ϵ8 = convert(T2, BigInt(-5097456320246635) // BigInt(161384393276391830))
    ϵ9 = convert(T2, BigInt(-14836286301942601) // BigInt(2200863462978094100))
    ϵ10 = convert(T2, BigInt(8712584401648310) // BigInt(230552839999042310))
    ϵ11 = convert(T2, BigInt(201598787147503) // BigInt(162515525987854870))
    ϵ12 = convert(T2, BigInt(-775224267694528) // BigInt(259632008240819930))
    ϵ13 = convert(T2, BigInt(37505781412258941025) // BigInt(343571135427505410))
    ϵ14 = convert(T2, BigInt(-70601687975788697909) // BigInt(197054662754227830))
    ϵ15 = convert(T2, BigInt(451638954762572460757303) // BigInt(1812919627060011327850))
    w1 = convert(T2, BigInt(922536264951379) // BigInt(23804051189793426))
    w8 = convert(T2, BigInt(225220497265063051) // BigInt(38410941914692596))
    w9 = convert(T2, BigInt(36148972660529108) // BigInt(25187799949529035))
    w10 = convert(T2, BigInt(-198426619086731179) // BigInt(28886469242593278))
    w11 = convert(T2, BigInt(6156310435334661) // BigInt(21948640431420694))
    w12 = convert(T2, BigInt(4487700950738411) // BigInt(30126758819889407))
    w13 = convert(T2, BigInt(29127285852344652967) // BigInt(26541231024453067))
    w14 = convert(T2, BigInt(-98614101020366183246) // BigInt(27482176375033227))
    w15 = convert(T2, BigInt(126863537044611133977) // BigInt(50939017338016615))
    w16 = convert(T2, BigInt(1421) // BigInt(3078))
    b21 = convert(T2, BigInt(3) // BigInt(184))
    b31 = convert(T2, BigInt(-433678958387722) // BigInt(60461976373487901))
    b32 = convert(T2, BigInt(800463031060200) // BigInt(17465287845577487))
    b41 = convert(T2, BigInt(377246711516913) // BigInt(26022242251007185))
    b43 = convert(T2, BigInt(1003144597305563) // BigInt(23065429003484804))
    b51 = convert(T2, BigInt(5232938947227754278) // BigInt(42414207887323079))
    b53 = convert(T2, BigInt(-18576057673326955826) // BigInt(47337109241323507))
    b54 = convert(T2, BigInt(9722802843766070268) // BigInt(36006158955671251))
    b61 = convert(T2, BigInt(1274983010403228) // BigInt(59501240010408527))
    b64 = convert(T2, BigInt(1653987885575953) // BigInt(24968953451991825))
    b65 = convert(T2, BigInt(70447607619) // BigInt(36673234034080855))
    b71 = convert(T2, BigInt(904721379210022) // BigInt(5952931408395025))
    b74 = convert(T2, BigInt(-23381308838657835) // BigInt(41704208630353019))
    b75 = convert(T2, BigInt(2066826469901) // BigInt(11925118158788070))
    b76 = convert(T2, BigInt(6895820618639281) // BigInt(11173940517165586))
    b81 = convert(T2, BigInt(703) // BigInt(25623))
    b86 = convert(T2, BigInt(2453249100410123) // BigInt(19386166204217813))
    b87 = convert(T2, BigInt(3412083729453488) // BigInt(36711207832699739))
    b91 = convert(T2, BigInt(13357) // BigInt(485888))
    b96 = convert(T2, BigInt(2029374314414594) // BigInt(16090898152406377))
    b97 = convert(T2, BigInt(1287982913423758) // BigInt(28873883427687143))
    b98 = convert(T2, BigInt(-6327) // BigInt(485888))
    b10_1 = convert(T2, BigInt(965192994961021) // BigInt(35149230657258373))
    b10_6 = convert(T2, BigInt(4127541420931126) // BigInt(32665647335015649))
    b10_7 = convert(T2, BigInt(3583837594502223) // BigInt(41664019148508749))
    b10_8 = convert(T2, BigInt(-246496779813877) // BigInt(27175788804018514))
    b10_9 = convert(T2, BigInt(28747455070668) // BigInt(8549877477842473))
    b11_1 = convert(T2, BigInt(63799716039495390) // BigInt(22202823855545599))
    b11_6 = convert(T2, BigInt(-43523) // BigInt(2744))
    b11_7 = convert(T2, BigInt(46511948996933887) // BigInt(19326005113913521))
    b11_8 = convert(T2, BigInt(4190224756003065953) // BigInt(25772283742746511))
    b11_9 = convert(T2, BigInt(4202436480827105093) // BigInt(56735195078669575))
    b11_10 = convert(T2, BigInt(-6231363253366161136) // BigInt(27638382159529459))
    b12_1 = convert(T2, BigInt(-2708709512927940344) // BigInt(43559765744627999))
    b12_6 = convert(T2, BigInt(9732832383152222769) // BigInt(28307424515467810))
    b12_7 = convert(T2, BigInt(-1467124908252536684) // BigInt(25216366962527245))
    b12_8 = convert(T2, BigInt(-123168960431653229617) // BigInt(37767835273866068))
    b12_9 = convert(T2, BigInt(-45716261318046911311) // BigInt(29725039922507753))
    b12_10 = convert(T2, BigInt(254588505817700334157) // BigInt(55659194024546537))
    b12_11 = convert(T2, BigInt(154479208885061992) // BigInt(61881010299032267))
    b13_1 = convert(T2, BigInt(7815889611024124839) // BigInt(10273038606706643))
    b13_6 = convert(T2, BigInt(-132865466831736468380) // BigInt(31621404109692561))
    b13_7 = convert(T2, BigInt(17703281196398756191) // BigInt(24701300463848133))
    b13_8 = convert(T2, BigInt(566474708394637174745) // BigInt(14201552570678562))
    b13_9 = convert(T2, BigInt(307285346872690815991) // BigInt(16345778788818343))
    b13_10 = convert(T2, BigInt(-2053676792180155828764) // BigInt(36714291312129259))
    b13_11 = convert(T2, BigInt(-78247540885035338) // BigInt(2948693499382457))
    b13_12 = convert(T2, BigInt(18260998947401049) // BigInt(15112168351518079))
    b14_1 = convert(T2, BigInt(21763183378499925071) // BigInt(28449009966588093))
    b14_6 = convert(T2, BigInt(-131859313072947008828) // BigInt(31210661767011655))
    b14_7 = convert(T2, BigInt(26844746596644780532) // BigInt(37252237041563455))
    b14_8 = convert(T2, BigInt(1253469865929367878849) // BigInt(31252862654010364))
    b14_9 = convert(T2, BigInt(463799691295259099617) // BigInt(24536685538656591))
    b14_10 = convert(T2, BigInt(-4389603511184674820033) // BigInt(78045747665396128))
    b14_11 = convert(T2, BigInt(-1128181715104199933) // BigInt(42276540809151831))
    b14_12 = convert(T2, BigInt(25502481921637862) // BigInt(20992853044100563))
    b14_13 = convert(T2, BigInt(-9559594005944) // BigInt(54330273344414971))
    b15_1 = convert(T2, BigInt(14646506801530214193) // BigInt(19100557410937333))
    b15_6 = convert(T2, BigInt(-100013811078882007482) // BigInt(23616683892616855))
    b15_7 = convert(T2, BigInt(51418842743714417428) // BigInt(71184230362278145))
    b15_8 = convert(T2, BigInt(449729579268755668521) // BigInt(11186469186029341))
    b15_9 = convert(T2, BigInt(830367055814329487149) // BigInt(43824976316629463))
    b15_10 = convert(T2, BigInt(-1510129653054461421096) // BigInt(26785760772223265))
    b15_11 = convert(T2, BigInt(-1401877877271794513) // BigInt(52404745175389990))
    b15_12 = convert(T2, BigInt(16736569978314491) // BigInt(13745172576768074))
    b15_13 = convert(T2, BigInt(-962836086693) // BigInt(28355588058859954))
    b15_14 = convert(T2, BigInt(-9941175715160) // BigInt(45586581308747583))
    b16_1 = convert(T2, BigInt(5917471745174541109) // BigInt(8344408594020065))
    b16_6 = convert(T2, BigInt(-107320434306581771105) // BigInt(27401651799311177))
    b16_7 = convert(T2, BigInt(5093729026420265343) // BigInt(7626561482131418))
    b16_8 = convert(T2, BigInt(2060040620258320258974) // BigInt(55399901444516827))
    b16_9 = convert(T2, BigInt(707106302890477200174) // BigInt(40350224587566143))
    b16_10 = convert(T2, BigInt(-1164323447229089993686) // BigInt(22328563534132263))
    b16_11 = convert(T2, BigInt(-831410876325262488) // BigInt(33593954156280041))
    b16_12 = convert(T2, BigInt(11423107778999224) // BigInt(9908693206689583))
    b16_13 = convert(T2, BigInt(165211034625068) // BigInt(39884924026051713))
    b16_14 = convert(T2, BigInt(-38387832271169) // BigInt(18010889018302554))

    QPRK98Tableau(d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14,
        ϵ1, ϵ8, ϵ9, ϵ10, ϵ11, ϵ12, ϵ13, ϵ14, ϵ15,
        w1, w8, w9, w10, w11, w12, w13, w14, w15, w16,
        b21, b31, b32, b41, b43, b51, b53, b54, b61, b64, b65, b71, b74, b75,
        b76, b81, b86, b87, b91, b96, b97, b98, b10_1, b10_6, b10_7, b10_8,
        b10_9, b11_1, b11_6, b11_7, b11_8, b11_9, b11_10, b12_1, b12_6, b12_7,
        b12_8, b12_9, b12_10, b12_11, b13_1, b13_6, b13_7, b13_8, b13_9, b13_10,
        b13_11, b13_12, b14_1, b14_6, b14_7, b14_8, b14_9, b14_10, b14_11, b14_12,
        b14_13, b15_1, b15_6, b15_7, b15_8, b15_9, b15_10, b15_11, b15_12, b15_13,
        b15_14, b16_1, b16_6, b16_7, b16_8, b16_9, b16_10, b16_11, b16_12, b16_13,
        b16_14)
end
