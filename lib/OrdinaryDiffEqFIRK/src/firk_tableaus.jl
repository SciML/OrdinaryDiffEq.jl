struct RadauIIA3Tableau{T, T2}
    T11::T
    T12::T
    T21::T
    T22::T
    TI11::T
    TI12::T
    TI21::T
    TI22::T
    c1::T2
    c2::T2
    α::T
    β::T
    e1::T
    e2::T
end

function RadauIIA3Tableau(T, T2)
    T11 = T(0.10540925533894596)
    T12 = T(-0.29814239699997197)
    T21 = T(0.9486832980505138)
    T22 = T(0.0)
    TI11 = T(0.0)
    TI12 = T(1.0540925533894598)
    TI21 = T(-3.3541019662496843)
    TI22 = T(0.3726779962499649)
    c1 = T2(1 / 3)
    c2 = T2(1.0)
    α = T(2.0)
    β = T(-sqrt(2))
    e1 = T(1 / 4)
    e2 = T(-1 / 4)
    RadauIIA3Tableau{T, T2}(T11, T12, T21, T22,
        TI11, TI12, TI21, TI22,
        c1, c2, α, β, e1, e2)
end

struct RadauIIA5Tableau{T, T2}
    T11::T
    T12::T
    T13::T
    T21::T
    T22::T
    T23::T
    T31::T
    #T32::T 
    #T33::T
    TI11::T
    TI12::T
    TI13::T
    TI21::T
    TI22::T
    TI23::T
    TI31::T
    TI32::T
    TI33::T
    c1::T2
    c2::T2
    #c3::T2
    γ::T
    α::T
    β::T
    e1::T
    e2::T
    e3::T
end

# inv(T) * inv(A) * T = [γ  0  0
#                        0  α -β
#                        0  β  α]
function RadauIIA5Tableau(T, T2)
    T11 = convert(T, 9.1232394870892942792e-02)
    T12 = convert(T, -0.14125529502095420843e0)
    T13 = convert(T, -3.0029194105147424492e-02)
    T21 = convert(T, 0.24171793270710701896e0)
    T22 = convert(T, 0.20412935229379993199e0)
    T23 = convert(T, 0.38294211275726193779e0)
    T31 = convert(T, 0.96604818261509293619e0)
    TI11 = convert(T, 4.3255798900631553510e0)
    TI12 = convert(T, 0.33919925181580986954e0)
    TI13 = convert(T, 0.54177053993587487119e0)
    TI21 = convert(T, -4.1787185915519047273e0)
    TI22 = convert(T, -0.32768282076106238708e0)
    TI23 = convert(T, 0.47662355450055045196e0)
    TI31 = convert(T, -0.50287263494578687595e0)
    TI32 = convert(T, 2.5719269498556054292e0)
    TI33 = convert(T, -0.59603920482822492497e0)

    sqrt6 = sqrt(6)

    c1 = convert(T2, (4 - sqrt6) / 10)
    c2 = convert(T2, (4 + sqrt6) / 10)

    cbrt9 = cbrt(9)

    γ′ = convert(T, (6.0 + cbrt9 * (cbrt9 - 1)) / 30) # eigval of `A`
    α′ = convert(T, (12.0 - cbrt9 * (cbrt9 - 1)) / 60) # eigval of `A`
    β′ = convert(T, cbrt9 * (cbrt9 + 1) * sqrt(3) / 60) # eigval of `A`
    scale = α′^2 + β′^2
    γ = inv(γ′)  # eigval of `inv(A)`
    α = α′ / scale # eigval of `inv(A)`
    β = β′ / scale # eigval of `inv(A)`

    e1 = convert(T, -(13 + 7 * sqrt6) / 3)
    e2 = convert(T, (-13 + 7 * sqrt6) / 3)
    e3 = convert(T, -1 / 3)
    RadauIIA5Tableau{T, T2}(T11, T12, T13, T21, T22, T23, T31, #= T33 = 0 =#
        TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33,
        c1, c2, #= c3 = 1 =#
        γ, α, β,
        e1, e2, e3)
end
#=
function BigRadauIIA5Tableau(T, T2)
    T11 = convert(T, 0.091232394870892942791548135249436196118684699372210280712184363514099824021240149574725365814781580305065489937969163922775110463056339192206701819661425186)
    T12 = convert(T, -0.141255295020954208427990383807797309409263248498594798844289981408804297900674604638610419147468875667691398225003133444988034605081071965848437945842767211)
    T31 = convert(T, -0.0300291941051474244918611170890538666683842974606300802563717702200388818691214144173874588956764952224874407424115249418136547481236684478531215095064078994)
    T21 = convert(T, 0.241717932707107018957474779310148232884879540532595279746187345714229132659465207414913313803429072060469564350914390845001169448350326344874859416624577348)
    T22 = convert(T, 0.204129352293799931995990810298338174086540402523315938937516234649384944528706774788799548853122282827246947911905379230680096946800308693162079538975632443)
    T23 = convert(T, 0.382942112757261937795438233599873210357792575012007744255205163027042915338009760005422153613194350161760232119048691964499888989151661861236831969497483828)
    T31 = convert(T, 0.966048182615092936190567080794590794996748754810883844283183333914131408744555961195911605614405476210484499875001737558078500322423463946527349731087504518)
    T32 = convert(T, 1.0)
    T33 = convert(T, 0.0)
    TI11 = convert(T, 4.32557989006315535102435095295614882731995158490590784287320458848019483341979047442263696495019938973156007686663488090615420049217658854859024016717169837)
    TI12 = convert(T, 0.339199251815809869542824974053410987511771566126056902312311333553438988409693737874718833892037643701271502187763370262948704203562215007824701228014200056)
    TI13 = convert(T, 0.541770539935874871186523033492089631898841317849243944095021379289933921771713116368931784890546144473788347538203807242114936998948954098533375649163016612)
    TI21 = convert(T, -4.17871859155190472734646265851205623000038388214686525896709481539843195209360778128456932548583273459040707932166364293012713818843609182148794380267482041)
    TI22 = convert(T, -0.327682820761062387082533272429616234245791838308340887801415258608836530255609335712523838667242449344879454518796849992049787172023800373390124427898159896)
    TI23 = convert(T, 0.476623554500550451960069084091012497939942928625055897109833707684876604712862299049343675491204859381277636585708398915065951363736337328178192801074535132)
    TI31 = convert(T, -0.502872634945786875951247343139544292859248429570937886791036339034110181540695221500843782634464164585836226038438397328726973424362168221527501738985822875)
    TI32 = convert(T, 2.57192694985560542918678535360167505469448742842178326395573566888176471664393761903447163100353067504020263109067033226021288356347565113471227052083596358)
    TI33 = convert(T, -0.596039204828224924968821911099302403289857517521591823052174732952989090998130905722763344484798508456930766594977798579939415052669401095404149917833710127)
    γ = convert(T, 3.63783425274449573220841851357777579794593608687391153215117488565841871456727143375130115708511223004183651123208497057248238260532214672028700625775335843)
    α = convert(T, 2.68108287362775213389579074321111210102703195656304423392441255717079064271636428312434942145744388497908174438395751471375880869733892663985649687112332242)
    β = convert(T, 3.05043019924741056942637762478756790444070419917947659226291744751211727051786694870515117615266028855554735929171362769761399150862332538376382934625577549)
    c1 = convert(T2, 0.155051025721682190180271592529410860803405251934332987156730743274903962254268497346014056689535976518140539877338581087514113454016224265837421604876272084)
    c2 = convert(T2, 0.644948974278317809819728407470589139196594748065667012843269256725096037745731502653985943310464023481859460122661418912485886545983775734162578395123729143)
    RadauIIA5Tableau{T, T2}(T11, T12, T13, T21, T22, T23, T31,T32, T33, 
    TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33,
    c1, c2, #= c3 = 1 =#
    γ, α, β,
    e1, e2, e3)
end
=#
struct RadauIIA9Tableau{T, T2}
    T11::T
    T12::T
    T13::T
    T14::T
    T15::T
    T21::T
    T22::T
    T23::T
    T24::T
    T25::T
    T31::T
    T32::T
    T33::T
    T34::T
    T35::T
    T41::T
    T42::T
    T43::T
    T44::T
    T45::T
    T51::T
    #T52::T
    #T53::T
    #T54::T
    #T55::T
    TI11::T
    TI12::T
    TI13::T
    TI14::T
    TI15::T
    TI21::T
    TI22::T
    TI23::T
    TI24::T
    TI25::T
    TI31::T
    TI32::T
    TI33::T
    TI34::T
    TI35::T
    TI41::T
    TI42::T
    TI43::T
    TI44::T
    TI45::T
    TI51::T
    TI52::T
    TI53::T
    TI54::T
    TI55::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    γ::T
    α1::T
    β1::T
    α2::T
    β2::T
    e1::T
    e2::T
    e3::T
    e4::T
    e5::T
end

function RadauIIA9Tableau(T, T2)
    T11 = convert(T, -1.251758622050104589014e-2)
    T12 = convert(T, -1.024204781790882707009e-2)
    T13 = convert(T, 4.767387729029572386318e-2)
    T14 = convert(T, -1.147851525522951470794e-2)
    T15 = convert(T, -1.401985889287541028108e-2)
    T21 = convert(T, -1.491670151895382429004e-3)
    T22 = convert(T, 5.017286451737105816299e-2)
    T23 = convert(T, -9.433181918161143698066e-2)
    T24 = convert(T, -7.668830749180162885157e-3)
    T25 = convert(T, 2.470857842651852681253e-2)
    T31 = convert(T, -7.298187638808714862266e-2)
    T32 = convert(T, -2.305395340434179467214e-1)
    T33 = convert(T, 1.027030453801258997922e-1)
    T34 = convert(T, 1.939846399882895091122e-2)
    T35 = convert(T, 8.180035370375117083639e-2)
    T41 = convert(T, -3.800914400035681041264e-1)
    T42 = convert(T, 3.778939022488612495439e-1)
    T43 = convert(T, 4.667441303324943592896e-1)
    T44 = convert(T, 4.076011712801990666217e-1)
    T45 = convert(T, 1.996824278868025259365e-1)
    T51 = convert(T, -9.219789736812104884883e-1)
    #T52 = convert(T, 1.000000000000000000000e0)
    #T53 = convert(T, 0.0000000000000000000000e0)
    #T54 = convert(T, 1.000000000000000000000e0)
    #T55 = convert(T, 0.0000000000000000000000e0)
    TI11 = convert(T, -3.004156772154440162771e1)
    TI12 = convert(T, -1.386510785627141316518e1)
    TI13 = convert(T, -3.480002774795185561828e0)
    TI14 = convert(T, 1.032008797825263422771e0)
    TI15 = convert(T, -8.043030450739899174753e-1)
    TI21 = convert(T, 5.344186437834911598895e0)
    TI22 = convert(T, 4.593615567759161004454e0)
    TI23 = convert(T, -3.036360323459424298646e0)
    TI24 = convert(T, 1.050660190231458863860e0)
    TI25 = convert(T, -2.727786118642962705386e-1)
    TI31 = convert(T, 3.748059807439804860051e0)
    TI32 = convert(T, -3.984965736343884667252e0)
    TI33 = convert(T, -1.044415641608018792942e0)
    TI34 = convert(T, 1.184098568137948487231e0)
    TI35 = convert(T, -4.499177701567803688988e-1)
    TI41 = convert(T, -3.304188021351900000806e1)
    TI42 = convert(T, -1.737695347906356701945e1)
    TI43 = convert(T, -1.721290632540055611515e-1)
    TI44 = convert(T, -9.916977798254264258817e-2)
    TI45 = convert(T, 5.312281158383066671849e-1)
    TI51 = convert(T, -8.611443979875291977700e0)
    TI52 = convert(T, 9.699991409528808231336e0)
    TI53 = convert(T, 1.914728639696874284851e0)
    TI54 = convert(T, 2.418692006084940026427e0)
    TI55 = convert(T, -1.047463487935337418694e0)

    c1 = convert(T2, 5.710419611451768219312e-2)
    c2 = convert(T2, 2.768430136381238276800e-1)
    c3 = convert(T2, 5.835904323689168200567e-1)
    c4 = convert(T2, 8.602401356562194478479e-1)
    #= c5 = convert(T2, 1) =#

    γ = convert(T, 6.286704751729276645173e0)
    α1 = convert(T, 3.655694325463572258243e0)
    β1 = convert(T, 6.543736899360077294021e0)
    α2 = convert(T, 5.700953298671789419170e0)
    β2 = convert(T, 3.210265600308549888425e0)

    e1 = convert(T, -2.778093394406463730479e1)
    e2 = convert(T, 3.641478498049213152712e0)
    e3 = convert(T, -1.252547721169118720491e0)
    e4 = convert(T, 5.920031671845428725662e-1)
    e5 = convert(T, -2.000000000000000000000e-1)

    RadauIIA9Tableau{T, T2}(T11, T12, T13, T14, T15,
        T21, T22, T23, T24, T25, T31, T32, T33, T34, T35,
        T41, T42, T43, T44, T45, T51, #=T52, T53, T54, T55=#
        TI11, TI12, TI13, TI14, TI15, TI21, TI22, TI23, TI24, TI25,
        TI31, TI32, TI33, TI34, TI35, TI41, TI42, TI43, TI44, TI45,
        TI51, TI52, TI53, TI54, TI55,
        c1, c2, c3, c4, #= c5 = 1 =#
        γ, α1, β1, α2, β2,
        e1, e2, e3, e4, e5)
end

struct adaptiveRadauTableau{T, T2, Int}
    T:: AbstractMatrix{T}
    TI::AbstractMatrix{T}
    γ::T
    α::AbstractVector{T}
    β::AbstractVector{T}
    c::AbstractVector{T}
    #e::AbstractVector{T}
    num_stages::Int
end

using Polynomials, GenericLinearAlgebra, LinearAlgebra, LinearSolve, GenericSchur, RootedTrees, Symbolics
using Symbolics: variables, variable, unwrap

function adaptiveRadauTableau(T, T2, num_stages::Int)
    tmp = Vector{BigFloat}(undef, num_stages - 1)
    for i in 1:(num_stages - 1)
        tmp[i] = 0
    end
    tmp2 = Vector{BigFloat}(undef, num_stages + 1)
    for i in 1:(num_stages + 1)
       tmp2[i]=(-1)^(num_stages + 1 - i) * binomial(num_stages , num_stages + 1 - i)
    end
    radau_p = Polynomial{BigFloat}([tmp; tmp2])
    for i in 1:(num_stages - 1)
        radau_p = derivative(radau_p)
    end
    c = roots(radau_p)
    c[num_stages] = 1
    c_powers = Matrix{BigFloat}(undef, num_stages, num_stages)
    for i in 1 : num_stages
        for j in 1 : num_stages
            c_powers[i,j] = c[i]^(j - 1)
        end
    end
    inverse_c_powers = inv(c_powers)
    c_q = Matrix{BigFloat}(undef, num_stages, num_stages)
    for i in 1 : num_stages
        for j in 1 : num_stages
            c_q[i,j] = c[i]^(j) / j
        end
    end
    a = c_q * inverse_c_powers
    a_inverse = inv(a)
    b = Vector{BigFloat}(undef, num_stages)
    for i in 1 : num_stages
        b[i] = a[num_stages, i]
    end
    vals = eigvals(a_inverse)
    γ = real(vals[num_stages])
    α = Vector{BigFloat}(undef, floor(Int, num_stages/2))
    β = Vector{BigFloat}(undef, floor(Int, num_stages/2))
    index = 1
    i = 1
    while i <= (num_stages - 1)
        α[index] = real(vals[i])
        β[index] = imag(vals[i + 1])
        index = index + 1
        i = i + 2
    end
    eigvec = eigvecs(a)
    vecs = Vector{Vector{BigFloat}}(undef, num_stages)
    i = 1
    index = 2
    while i < num_stages 
        vecs[index] = real(eigvec[:, i] ./ eigvec[num_stages, i])
        vecs[index + 1] = -imag(eigvec[:, i] ./ eigvec[num_stages, i])
        index += 2
        i += 2
    end
    vecs[1] = real(eigvec[:, num_stages])
    tmp3 = vcat(vecs)
    T = Matrix{BigFloat}(undef, num_stages, num_stages)
    for j in 1 : num_stages
        for i in 1 : num_stages
            T[i, j] = tmp3[j][i]
        end
    end
    TI = inv(T)
    #=
    p = num_stages
    eb = variables(:b, 1:num_stages + 1)
    @variables y
    zz = zeros(size(a, 1) + 1)
    zz2 = zeros(size(a, 1))
    eA = [zz'
          zz2 a]
    ec = [0; c]
    constraints = map(Iterators.flatten(RootedTreeIterator(i) for i in 1:p)) do t
        residual_order_condition(t, RungeKuttaMethod(eA, eb, ec))
    end
    AA, bb, islinear = Symbolics.linear_expansion(substitute.(constraints, (eb[1]=>y,)), eb[2:end])
    AA = BigFloat.(map(unwrap, AA))
    idxs = qr(AA', ColumnNorm()).p[1:num_stages]
    @assert rank(AA[idxs, :]) == num_stages
    @assert islinear
    Symbolics.expand.((AA[idxs, :] \ -bb[idxs]) - b)=#
    #e = b_hat - b
    adaptiveRadauTableau{Any, T2, Int}(T, TI, γ, α, β, c, num_stages)
end