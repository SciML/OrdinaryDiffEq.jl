@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern6Cache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 9 || always_calc_begin
        (; c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98) = cache.tab
        (; k1, k2, k3, k4, k5, k6, k7, k8, k9, tmp) = cache
        @.. broadcast = false tmp = uprev + dt * (a21 * k1)
        f(k2, tmp, p, t + c1 * dt)
        @.. broadcast = false tmp = uprev + dt * (a31 * k1 + a32 * k2)
        f(k3, tmp, p, t + c2 * dt)
        @.. broadcast = false tmp = uprev + dt * (a41 * k1 + a43 * k3)
        f(k4, tmp, p, t + c3 * dt)
        @.. broadcast = false tmp = uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
        f(k5, tmp, p, t + c4 * dt)
        @.. broadcast = false tmp = uprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
        f(k6, tmp, p, t + c5 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
        f(k7, tmp, p, t + c6 * dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 +
                a87 * k7
        )
        f(k8, tmp, p, t + dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a91 * k1 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 +
                a98 * k8
        )
        f(k9, tmp, p, t + dt)
        copyat_or_push!(k, 1, k1)
        copyat_or_push!(k, 2, k2)
        copyat_or_push!(k, 3, k3)
        copyat_or_push!(k, 4, k4)
        copyat_or_push!(k, 5, k5)
        copyat_or_push!(k, 6, k6)
        copyat_or_push!(k, 7, k7)
        copyat_or_push!(k, 8, k8)
        copyat_or_push!(k, 9, k9)
    end
    if (allow_calc_end && length(k) < 12) || force_calc_end # Have not added the extra stages yet
        (; c10, a1001, a1004, a1005, a1006, a1007, a1008, a1009, c11, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, c12, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211) = cache.tab.extra
        (; tmp) = cache
        rtmp = similar(cache.k1)
        uidx = eachindex(uprev)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1001 * k[1] + a1004 * k[4] + a1005 * k[5] + a1006 * k[6] +
                a1007 * k[7] + a1008 * k[8] + a1009 * k[9]
        )
        f(rtmp, tmp, p, t + c10 * dt)
        copyat_or_push!(k, 10, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                a1107 * k[7] + a1108 * k[8] + a1109 * k[9] + a1110 * k[10]
        )
        f(rtmp, tmp, p, t + c11 * dt)
        copyat_or_push!(k, 11, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                a1207 * k[7] + a1208 * k[8] + a1209 * k[9] +
                a1210 * k[10] + a1211 * k[11]
        )
        f(rtmp, tmp, p, t + c12 * dt)
        copyat_or_push!(k, 12, rtmp)
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern7CacheType,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    if length(k) < 10 || always_calc_begin
        @OnDemandTableauExtract Vern7Tableau T T2
        (; k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, tmp) = cache
        f(k1, uprev, p, t)
        @.. broadcast = false tmp = uprev + dt * (a021 * k1)
        f(k2, tmp, p, t + c2 * dt)
        @.. broadcast = false tmp = uprev + dt * (a031 * k1 + a032 * k2)
        f(k3, tmp, p, t + c3 * dt)
        @.. broadcast = false tmp = uprev + dt * (a041 * k1 + a043 * k3)
        f(k4, tmp, p, t + c4 * dt)
        @.. broadcast = false tmp = uprev + dt * (a051 * k1 + a053 * k3 + a054 * k4)
        f(k5, tmp, p, t + c5 * dt)
        @.. broadcast = false tmp = uprev + dt * (a061 * k1 + a063 * k3 + a064 * k4 + a065 * k5)
        f(k6, tmp, p, t + c6 * dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (a071 * k1 + a073 * k3 + a074 * k4 + a075 * k5 + a076 * k6)
        f(k7, tmp, p, t + c7 * dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a081 * k1 + a083 * k3 + a084 * k4 + a085 * k5 + a086 * k6 +
                a087 * k7
        )
        f(k8, tmp, p, t + c8 * dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a091 * k1 + a093 * k3 + a094 * k4 + a095 * k5 + a096 * k6 +
                a097 * k7 + a098 * k8
        )
        f(k9, tmp, p, t + dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a101 * k1 + a103 * k3 + a104 * k4 + a105 * k5 + a106 * k6 +
                a107 * k7
        )
        f(k10, tmp, p, t + dt)
        copyat_or_push!(k, 1, k1)
        copyat_or_push!(k, 2, k2)
        copyat_or_push!(k, 3, k3)
        copyat_or_push!(k, 4, k4)
        copyat_or_push!(k, 5, k5)
        copyat_or_push!(k, 6, k6)
        copyat_or_push!(k, 7, k7)
        copyat_or_push!(k, 8, k8)
        copyat_or_push!(k, 9, k9)
        copyat_or_push!(k, 10, k10)
    end
    if (allow_calc_end && length(k) < 16) || force_calc_end # Have not added the extra stages yet
        (; tmp) = cache
        rtmp = similar(cache.k1)
        @OnDemandTableauExtract Vern7ExtraStages T T2
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                a1107 * k[7] + a1108 * k[8] + a1109 * k[9]
        )
        f(rtmp, tmp, p, t + c11 * dt)
        copyat_or_push!(k, 11, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1211 * k[11]
        )
        f(rtmp, tmp, p, t + c12 * dt)
        copyat_or_push!(k, 12, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1301 * k[1] + a1304 * k[4] + a1305 * k[5] + a1306 * k[6] +
                a1307 * k[7] + a1308 * k[8] + a1309 * k[9] +
                a1311 * k[11] + a1312 * k[12]
        )
        f(rtmp, tmp, p, t + c13 * dt)
        copyat_or_push!(k, 13, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1401 * k[1] + a1404 * k[4] + a1405 * k[5] + a1406 * k[6] +
                a1407 * k[7] + a1408 * k[8] + a1409 * k[9] +
                a1411 * k[11] + a1412 * k[12] + a1413 * k[13]
        )
        f(rtmp, tmp, p, t + c14 * dt)
        copyat_or_push!(k, 14, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1501 * k[1] + a1504 * k[4] + a1505 * k[5] + a1506 * k[6] +
                a1507 * k[7] + a1508 * k[8] + a1509 * k[9] +
                a1511 * k[11] + a1512 * k[12] + a1513 * k[13]
        )
        f(rtmp, tmp, p, t + c15 * dt)
        copyat_or_push!(k, 15, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1601 * k[1] + a1604 * k[4] + a1605 * k[5] + a1606 * k[6] +
                a1607 * k[7] + a1608 * k[8] + a1609 * k[9] +
                a1611 * k[11] + a1612 * k[12] + a1613 * k[13]
        )
        f(rtmp, tmp, p, t + c16 * dt)
        copyat_or_push!(k, 16, rtmp)
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern8Cache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 13 || always_calc_begin
        (; c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701, a0704, a0705, a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906, a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1304, a1305, a1306, a1307, a1308, a1309, a1310) = cache.tab
        (; k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, tmp) = cache
        f(k1, uprev, p, t)
        @.. broadcast = false tmp = uprev + dt * (a0201 * k1)
        f(k2, tmp, p, t + c2 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0301 * k1 + a0302 * k2)
        f(k3, tmp, p, t + c3 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0401 * k1 + a0403 * k3)
        f(k4, tmp, p, t + c4 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0501 * k1 + a0503 * k3 + a0504 * k4)
        f(k5, tmp, p, t + c5 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0601 * k1 + a0604 * k4 + a0605 * k5)
        f(k6, tmp, p, t + c6 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6)
        f(k7, tmp, p, t + c7 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a0801 * k1 + a0804 * k4 + a0805 * k5 + a0806 * k6 +
                a0807 * k7
        )
        f(k8, tmp, p, t + c8 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a0901 * k1 + a0904 * k4 + a0905 * k5 + a0906 * k6 +
                a0907 * k7 + a0908 * k8
        )
        f(k9, tmp, p, t + c9 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1001 * k1 + a1004 * k4 + a1005 * k5 + a1006 * k6 +
                a1007 * k7 + a1008 * k8 + a1009 * k9
        )
        f(k10, tmp, p, t + c10 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1101 * k1 + a1104 * k4 + a1105 * k5 + a1106 * k6 +
                a1107 * k7 + a1108 * k8 + a1109 * k9 + a1110 * k10
        )
        f(k11, tmp, p, t + c11 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1201 * k1 + a1204 * k4 + a1205 * k5 + a1206 * k6 +
                a1207 * k7 + a1208 * k8 + a1209 * k9 + a1210 * k10 +
                a1211 * k11
        )
        f(k12, tmp, p, t + dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1301 * k1 + a1304 * k4 + a1305 * k5 + a1306 * k6 +
                a1307 * k7 + a1308 * k8 + a1309 * k9 + a1310 * k10
        )
        f(k13, tmp, p, t + dt)
        copyat_or_push!(k, 1, k1)
        copyat_or_push!(k, 2, k2)
        copyat_or_push!(k, 3, k3)
        copyat_or_push!(k, 4, k4)
        copyat_or_push!(k, 5, k5)
        copyat_or_push!(k, 6, k6)
        copyat_or_push!(k, 7, k7)
        copyat_or_push!(k, 8, k8)
        copyat_or_push!(k, 9, k9)
        copyat_or_push!(k, 10, k10)
        copyat_or_push!(k, 11, k11)
        copyat_or_push!(k, 12, k12)
        copyat_or_push!(k, 13, k13)
    end
    if (allow_calc_end && length(k) < 21) || force_calc_end # Have not added the extra stages yet
        rtmp = similar(cache.k1)
        (; c14, a1401, a1406, a1407, a1408, a1409, a1410, a1411, a1412, c15, a1501, a1506, a1507, a1508, a1509, a1510, a1511, a1512, a1514, c16, a1601, a1606, a1607, a1608, a1609, a1610, a1611, a1612, a1614, a1615, c17, a1701, a1706, a1707, a1708, a1709, a1710, a1711, a1712, a1714, a1715, a1716, c18, a1801, a1806, a1807, a1808, a1809, a1810, a1811, a1812, a1814, a1815, a1816, a1817, c19, a1901, a1906, a1907, a1908, a1909, a1910, a1911, a1912, a1914, a1915, a1916, a1917, c20, a2001, a2006, a2007, a2008, a2009, a2010, a2011, a2012, a2014, a2015, a2016, a2017, c21, a2101, a2106, a2107, a2108, a2109, a2110, a2111, a2112, a2114, a2115, a2116, a2117) = cache.tab.extra
        (; tmp) = cache
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1401 * k[1] + a1406 * k[6] + a1407 * k[7] + a1408 * k[8] +
                a1409 * k[9] + a1410 * k[10] + a1411 * k[11] +
                a1412 * k[12]
        )
        f(rtmp, tmp, p, t + c14 * dt)
        copyat_or_push!(k, 14, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1501 * k[1] + a1506 * k[6] + a1507 * k[7] + a1508 * k[8] +
                a1509 * k[9] + a1510 * k[10] + a1511 * k[11] +
                a1512 * k[12] + a1514 * k[14]
        )
        f(rtmp, tmp, p, t + c15 * dt)
        copyat_or_push!(k, 15, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1601 * k[1] + a1606 * k[6] + a1607 * k[7] + a1608 * k[8] +
                a1609 * k[9] + a1610 * k[10] + a1611 * k[11] +
                a1612 * k[12] + a1614 * k[14] + a1615 * k[15]
        )
        f(rtmp, tmp, p, t + c16 * dt)
        copyat_or_push!(k, 16, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1701 * k[1] + a1706 * k[6] + a1707 * k[7] + a1708 * k[8] +
                a1709 * k[9] + a1710 * k[10] + a1711 * k[11] +
                a1712 * k[12] + a1714 * k[14] + a1715 * k[15] +
                a1716 * k[16]
        )
        f(rtmp, tmp, p, t + c17 * dt)
        copyat_or_push!(k, 17, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1801 * k[1] + a1806 * k[6] + a1807 * k[7] + a1808 * k[8] +
                a1809 * k[9] + a1810 * k[10] + a1811 * k[11] +
                a1812 * k[12] + a1814 * k[14] + a1815 * k[15] +
                a1816 * k[16] + a1817 * k[17]
        )
        f(rtmp, tmp, p, t + c18 * dt)
        copyat_or_push!(k, 18, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1901 * k[1] + a1906 * k[6] + a1907 * k[7] + a1908 * k[8] +
                a1909 * k[9] + a1910 * k[10] + a1911 * k[11] +
                a1912 * k[12] + a1914 * k[14] + a1915 * k[15] +
                a1916 * k[16] + a1917 * k[17]
        )
        f(rtmp, tmp, p, t + c19 * dt)
        copyat_or_push!(k, 19, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2001 * k[1] + a2006 * k[6] + a2007 * k[7] + a2008 * k[8] +
                a2009 * k[9] + a2010 * k[10] + a2011 * k[11] +
                a2012 * k[12] + a2014 * k[14] + a2015 * k[15] +
                a2016 * k[16] + a2017 * k[17]
        )
        f(rtmp, tmp, p, t + c20 * dt)
        copyat_or_push!(k, 20, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2101 * k[1] + a2106 * k[6] + a2107 * k[7] + a2108 * k[8] +
                a2109 * k[9] + a2110 * k[10] + a2111 * k[11] +
                a2112 * k[12] + a2114 * k[14] + a2115 * k[15] +
                a2116 * k[16] + a2117 * k[17]
        )
        f(rtmp, tmp, p, t + c21 * dt)
        copyat_or_push!(k, 21, rtmp)
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern9Cache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    if length(k) < 10 || always_calc_begin
        @OnDemandTableauExtract Vern9Tableau T T2
        (; k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, tmp) = cache
        uidx = eachindex(uprev)
        f(k1, uprev, p, t)
        @.. broadcast = false tmp = uprev + dt * (a0201 * k1)
        f(k2, tmp, p, t + c1 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0301 * k1 + a0302 * k2)
        f(k3, tmp, p, t + c2 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0401 * k1 + a0403 * k3)
        f(k4, tmp, p, t + c3 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0501 * k1 + a0503 * k3 + a0504 * k4)
        f(k5, tmp, p, t + c4 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0601 * k1 + a0604 * k4 + a0605 * k5)
        f(k6, tmp, p, t + c5 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6)
        f(k7, tmp, p, t + c6 * dt)
        @.. broadcast = false tmp = uprev + dt * (a0801 * k1 + a0806 * k6 + a0807 * k7)
        f(k8, tmp, p, t + c7 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (a0901 * k1 + a0906 * k6 + a0907 * k7 + a0908 * k8)
        f(k9, tmp, p, t + c8 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1001 * k1 + a1006 * k6 + a1007 * k7 + a1008 * k8 +
                a1009 * k9
        )
        f(k10, tmp, p, t + c9 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1101 * k1 + a1106 * k6 + a1107 * k7 + a1108 * k8 +
                a1109 * k9 + a1110 * k10
        )
        f(k11, tmp, p, t + c10 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1201 * k1 + a1206 * k6 + a1207 * k7 + a1208 * k8 +
                a1209 * k9 + a1210 * k10 + a1211 * k11
        )
        f(k12, tmp, p, t + c11 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1301 * k1 + a1306 * k6 + a1307 * k7 + a1308 * k8 +
                a1309 * k9 + a1310 * k10 + a1311 * k11 + a1312 * k12
        )
        f(k13, tmp, p, t + c12 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1401 * k1 + a1406 * k6 + a1407 * k7 + a1408 * k8 +
                a1409 * k9 + a1410 * k10 + a1411 * k11 + a1412 * k12 +
                a1413 * k13
        )
        f(k14, tmp, p, t + c13 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1501 * k1 + a1506 * k6 + a1507 * k7 + a1508 * k8 +
                a1509 * k9 + a1510 * k10 + a1511 * k11 + a1512 * k12 +
                a1513 * k13 + a1514 * k14
        )
        f(k15, tmp, p, t + dt)
        @.. broadcast = false tmp = uprev +
            dt * (
            a1601 * k1 + a1606 * k6 + a1607 * k7 + a1608 * k8 +
                a1609 * k9 + a1610 * k10 + a1611 * k11 + a1612 * k12 +
                a1613 * k13
        )
        f(k16, tmp, p, t + dt)
        copyat_or_push!(k, 1, k1)
        copyat_or_push!(k, 2, k8)
        copyat_or_push!(k, 3, k9)
        copyat_or_push!(k, 4, k10)
        copyat_or_push!(k, 5, k11)
        copyat_or_push!(k, 6, k12)
        copyat_or_push!(k, 7, k13)
        copyat_or_push!(k, 8, k14)
        copyat_or_push!(k, 9, k15)
        copyat_or_push!(k, 10, k16)
    end
    if (allow_calc_end && length(k) < 20) || force_calc_end # Have not added the extra stages yet
        rtmp = similar(cache.k1)
        uidx = eachindex(uprev)
        (; tmp) = cache
        @OnDemandTableauExtract Vern9ExtraStages T T2
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1701 * k[1] + a1708 * k[2] + a1709 * k[3] + a1710 * k[4] +
                a1711 * k[5] + a1712 * k[6] + a1713 * k[7] + a1714 * k[8] +
                a1715 * k[9]
        )
        f(rtmp, tmp, p, t + c17 * dt)
        copyat_or_push!(k, 11, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1801 * k[1] + a1808 * k[2] + a1809 * k[3] + a1810 * k[4] +
                a1811 * k[5] + a1812 * k[6] + a1813 * k[7] + a1814 * k[8] +
                a1815 * k[9] + a1817 * k[11]
        )
        f(rtmp, tmp, p, t + c18 * dt)
        copyat_or_push!(k, 12, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a1901 * k[1] + a1908 * k[2] + a1909 * k[3] + a1910 * k[4] +
                a1911 * k[5] + a1912 * k[6] + a1913 * k[7] + a1914 * k[8] +
                a1915 * k[9] + a1917 * k[11] + a1918 * k[12]
        )
        f(rtmp, tmp, p, t + c19 * dt)
        copyat_or_push!(k, 13, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2001 * k[1] + a2008 * k[2] + a2009 * k[3] + a2010 * k[4] +
                a2011 * k[5] + a2012 * k[6] + a2013 * k[7] + a2014 * k[8] +
                a2015 * k[9] + a2017 * k[11] + a2018 * k[12] +
                a2019 * k[13]
        )
        f(rtmp, tmp, p, t + c20 * dt)
        copyat_or_push!(k, 14, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2101 * k[1] + a2108 * k[2] + a2109 * k[3] + a2110 * k[4] +
                a2111 * k[5] + a2112 * k[6] + a2113 * k[7] + a2114 * k[8] +
                a2115 * k[9] + a2117 * k[11] + a2118 * k[12] +
                a2119 * k[13] + a2120 * k[14]
        )
        f(rtmp, tmp, p, t + c21 * dt)
        copyat_or_push!(k, 15, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2201 * k[1] + a2208 * k[2] + a2209 * k[3] + a2210 * k[4] +
                a2211 * k[5] + a2212 * k[6] + a2213 * k[7] + a2214 * k[8] +
                a2215 * k[9] + a2217 * k[11] + a2218 * k[12] +
                a2219 * k[13] + a2220 * k[14] + a2221 * k[15]
        )
        f(rtmp, tmp, p, t + c22 * dt)
        copyat_or_push!(k, 16, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2301 * k[1] + a2308 * k[2] + a2309 * k[3] + a2310 * k[4] +
                a2311 * k[5] + a2312 * k[6] + a2313 * k[7] + a2314 * k[8] +
                a2315 * k[9] + a2317 * k[11] + a2318 * k[12] +
                a2319 * k[13] + a2320 * k[14] + a2321 * k[15]
        )
        f(rtmp, tmp, p, t + c23 * dt)
        copyat_or_push!(k, 17, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2401 * k[1] + a2408 * k[2] + a2409 * k[3] + a2410 * k[4] +
                a2411 * k[5] + a2412 * k[6] + a2413 * k[7] + a2414 * k[8] +
                a2415 * k[9] + a2417 * k[11] + a2418 * k[12] +
                a2419 * k[13] + a2420 * k[14] + a2421 * k[15]
        )
        f(rtmp, tmp, p, t + c24 * dt)
        copyat_or_push!(k, 18, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2501 * k[1] + a2508 * k[2] + a2509 * k[3] + a2510 * k[4] +
                a2511 * k[5] + a2512 * k[6] + a2513 * k[7] + a2514 * k[8] +
                a2515 * k[9] + a2517 * k[11] + a2518 * k[12] +
                a2519 * k[13] + a2520 * k[14] + a2521 * k[15]
        )
        f(rtmp, tmp, p, t + c25 * dt)
        copyat_or_push!(k, 19, rtmp)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a2601 * k[1] + a2608 * k[2] + a2609 * k[3] + a2610 * k[4] +
                a2611 * k[5] + a2612 * k[6] + a2613 * k[7] + a2614 * k[8] +
                a2615 * k[9] + a2617 * k[11] + a2618 * k[12] +
                a2619 * k[13] + a2620 * k[14] + a2621 * k[15]
        )
        f(rtmp, tmp, p, t + c26 * dt)
        copyat_or_push!(k, 20, rtmp)
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern6ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 9 || always_calc_begin
        (; c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98) = cache.tab
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(uprev + dt * (a21 * k[1]), p, t + c1 * dt))
        copyat_or_push!(k, 3, f(uprev + dt * (a31 * k[1] + a32 * k[2]), p, t + c2 * dt))
        copyat_or_push!(k, 4, f(uprev + dt * (a41 * k[1] + a43 * k[3]), p, t + c3 * dt))
        copyat_or_push!(
            k, 5,
            f(
                uprev + dt * (a51 * k[1] + a53 * k[3] + a54 * k[4]), p,
                t + c4 * dt
            )
        )
        copyat_or_push!(
            k, 6,
            f(
                uprev + dt * (a61 * k[1] + a63 * k[3] + a64 * k[4] + a65 * k[5]),
                p, t + c5 * dt
            )
        )
        copyat_or_push!(
            k, 7,
            f(
                uprev +
                    dt *
                    (a71 * k[1] + a73 * k[3] + a74 * k[4] + a75 * k[5] + a76 * k[6]),
                p, t + c6 * dt
            )
        )
        copyat_or_push!(
            k, 8,
            f(
                uprev +
                    dt *
                    (
                    a81 * k[1] + a83 * k[3] + a84 * k[4] + a85 * k[5] + a86 * k[6] +
                        a87 * k[7]
                ),
                p,
                t + dt
            )
        )
        copyat_or_push!(
            k, 9,
            f(
                uprev +
                    dt *
                    (
                    a91 * k[1] + a94 * k[4] + a95 * k[5] + a96 * k[6] + a97 * k[7] +
                        a98 * k[8]
                ),
                p,
                t + dt
            )
        )
    end
    if (allow_calc_end && length(k) < 12) || force_calc_end # Have not added the extra stages yet
        (; c10, a1001, a1004, a1005, a1006, a1007, a1008, a1009, c11, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, c12, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211) = cache.tab.extra
        copyat_or_push!(
            k, 10,
            f(
                uprev +
                    dt *
                    (
                    a1001 * k[1] + a1004 * k[4] + a1005 * k[5] + a1006 * k[6] +
                        a1007 * k[7] + a1008 * k[8] + a1009 * k[9]
                ),
                p,
                t + c10 * dt
            )
        )
        copyat_or_push!(
            k, 11,
            f(
                uprev +
                    dt *
                    (
                    a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                        a1107 * k[7] + a1108 * k[8] + a1109 * k[9] + a1110 * k[10]
                ),
                p,
                t + c11 * dt
            )
        )
        copyat_or_push!(
            k, 12,
            f(
                uprev +
                    dt *
                    (
                    a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                        a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1210 * k[10] +
                        a1211 * k[11]
                ),
                p,
                t + c12 * dt
            )
        )
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern7ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    if length(k) < 10 || always_calc_begin
        @OnDemandTableauExtract Vern7Tableau T T2
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(uprev + dt * (a021 * k[1]), p, t + c2 * dt))
        copyat_or_push!(k, 3, f(uprev + dt * (a031 * k[1] + a032 * k[2]), p, t + c3 * dt))
        copyat_or_push!(k, 4, f(uprev + dt * (a041 * k[1] + a043 * k[3]), p, t + c4 * dt))
        copyat_or_push!(
            k, 5,
            f(
                uprev + dt * (a051 * k[1] + a053 * k[3] + a054 * k[4]), p,
                t + c5 * dt
            )
        )
        copyat_or_push!(
            k, 6,
            f(
                uprev +
                    dt * (a061 * k[1] + a063 * k[3] + a064 * k[4] + a065 * k[5]), p,
                t + c6 * dt
            )
        )
        copyat_or_push!(
            k, 7,
            f(
                uprev +
                    dt * (
                    a071 * k[1] + a073 * k[3] + a074 * k[4] + a075 * k[5] +
                        a076 * k[6]
                ),
                p,
                t + c7 * dt
            )
        )
        copyat_or_push!(
            k, 8,
            f(
                uprev +
                    dt * (
                    a081 * k[1] + a083 * k[3] + a084 * k[4] + a085 * k[5] +
                        a086 * k[6] + a087 * k[7]
                ),
                p,
                t + c8 * dt
            )
        )
        copyat_or_push!(
            k, 9,
            f(
                uprev +
                    dt * (
                    a091 * k[1] + a093 * k[3] + a094 * k[4] + a095 * k[5] +
                        a096 * k[6] + a097 * k[7] + a098 * k[8]
                ),
                p,
                t + dt
            )
        )
        copyat_or_push!(
            k, 10,
            f(
                uprev +
                    dt * (
                    a101 * k[1] + a103 * k[3] + a104 * k[4] + a105 * k[5] +
                        a106 * k[6] + a107 * k[7]
                ),
                p,
                t + dt
            )
        )
    end
    if (allow_calc_end && length(k) < 16) || force_calc_end # Have not added the extra stages yet
        @OnDemandTableauExtract Vern7ExtraStages T T2
        copyat_or_push!(
            k, 11,
            f(
                uprev +
                    dt *
                    (
                    a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                        a1107 * k[7] + a1108 * k[8] + a1109 * k[9]
                ),
                p,
                t + c11 * dt
            )
        )
        copyat_or_push!(
            k, 12,
            f(
                uprev +
                    dt *
                    (
                    a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                        a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1211 * k[11]
                ),
                p,
                t + c12 * dt
            )
        )
        copyat_or_push!(
            k, 13,
            f(
                uprev +
                    dt *
                    (
                    a1301 * k[1] + a1304 * k[4] + a1305 * k[5] + a1306 * k[6] +
                        a1307 * k[7] + a1308 * k[8] + a1309 * k[9] + a1311 * k[11] +
                        a1312 * k[12]
                ),
                p,
                t + c13 * dt
            )
        )
        copyat_or_push!(
            k, 14,
            f(
                uprev +
                    dt *
                    (
                    a1401 * k[1] + a1404 * k[4] + a1405 * k[5] + a1406 * k[6] +
                        a1407 * k[7] + a1408 * k[8] + a1409 * k[9] + a1411 * k[11] +
                        a1412 * k[12] + a1413 * k[13]
                ),
                p,
                t + c14 * dt
            )
        )
        copyat_or_push!(
            k, 15,
            f(
                uprev +
                    dt *
                    (
                    a1501 * k[1] + a1504 * k[4] + a1505 * k[5] + a1506 * k[6] +
                        a1507 * k[7] + a1508 * k[8] + a1509 * k[9] + a1511 * k[11] +
                        a1512 * k[12] + a1513 * k[13]
                ),
                p,
                t + c15 * dt
            )
        )
        copyat_or_push!(
            k, 16,
            f(
                uprev +
                    dt *
                    (
                    a1601 * k[1] + a1604 * k[4] + a1605 * k[5] + a1606 * k[6] +
                        a1607 * k[7] + a1608 * k[8] + a1609 * k[9] + a1611 * k[11] +
                        a1612 * k[12] + a1613 * k[13]
                ),
                p,
                t + c16 * dt
            )
        )
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern8ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 13 || always_calc_begin
        (; c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701, a0704, a0705, a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906, a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1304, a1305, a1306, a1307, a1308, a1309, a1310) = cache.tab
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(uprev + dt * (a0201 * k[1]), p, t + c2 * dt))
        copyat_or_push!(k, 3, f(uprev + dt * (a0301 * k[1] + a0302 * k[2]), p, t + c3 * dt))
        copyat_or_push!(k, 4, f(uprev + dt * (a0401 * k[1] + a0403 * k[3]), p, t + c4 * dt))
        copyat_or_push!(
            k, 5,
            f(
                uprev + dt * (a0501 * k[1] + a0503 * k[3] + a0504 * k[4]), p,
                t + c5 * dt
            )
        )
        copyat_or_push!(
            k, 6,
            f(
                uprev + dt * (a0601 * k[1] + a0604 * k[4] + a0605 * k[5]), p,
                t + c6 * dt
            )
        )
        copyat_or_push!(
            k, 7,
            f(
                uprev +
                    dt * (a0701 * k[1] + a0704 * k[4] + a0705 * k[5] + a0706 * k[6]),
                p, t + c7 * dt
            )
        )
        copyat_or_push!(
            k, 8,
            f(
                uprev +
                    dt *
                    (
                    a0801 * k[1] + a0804 * k[4] + a0805 * k[5] + a0806 * k[6] +
                        a0807 * k[7]
                ),
                p,
                t + c8 * dt
            )
        )
        copyat_or_push!(
            k, 9,
            f(
                uprev +
                    dt *
                    (
                    a0901 * k[1] + a0904 * k[4] + a0905 * k[5] + a0906 * k[6] +
                        a0907 * k[7] + a0908 * k[8]
                ),
                p,
                t + c9 * dt
            )
        )
        copyat_or_push!(
            k, 10,
            f(
                uprev +
                    dt *
                    (
                    a1001 * k[1] + a1004 * k[4] + a1005 * k[5] + a1006 * k[6] +
                        a1007 * k[7] + a1008 * k[8] + a1009 * k[9]
                ),
                p,
                t + c10 * dt
            )
        )
        copyat_or_push!(
            k, 11,
            f(
                uprev +
                    dt *
                    (
                    a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                        a1107 * k[7] + a1108 * k[8] + a1109 * k[9] + a1110 * k[10]
                ),
                p,
                t + c11 * dt
            )
        )
        copyat_or_push!(
            k, 12,
            f(
                uprev +
                    dt *
                    (
                    a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                        a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1210 * k[10] +
                        a1211 * k[11]
                ),
                p,
                t + dt
            )
        )
        copyat_or_push!(
            k, 13,
            f(
                uprev +
                    dt *
                    (
                    a1301 * k[1] + a1304 * k[4] + a1305 * k[5] + a1306 * k[6] +
                        a1307 * k[7] + a1308 * k[8] + a1309 * k[9] + a1310 * k[10]
                ),
                p,
                t + dt
            )
        )
    end
    if (allow_calc_end && length(k) < 21) || force_calc_end # Have not added the extra stages yet
        (; c14, a1401, a1406, a1407, a1408, a1409, a1410, a1411, a1412, c15, a1501, a1506, a1507, a1508, a1509, a1510, a1511, a1512, a1514, c16, a1601, a1606, a1607, a1608, a1609, a1610, a1611, a1612, a1614, a1615, c17, a1701, a1706, a1707, a1708, a1709, a1710, a1711, a1712, a1714, a1715, a1716, c18, a1801, a1806, a1807, a1808, a1809, a1810, a1811, a1812, a1814, a1815, a1816, a1817, c19, a1901, a1906, a1907, a1908, a1909, a1910, a1911, a1912, a1914, a1915, a1916, a1917, c20, a2001, a2006, a2007, a2008, a2009, a2010, a2011, a2012, a2014, a2015, a2016, a2017, c21, a2101, a2106, a2107, a2108, a2109, a2110, a2111, a2112, a2114, a2115, a2116, a2117) = cache.tab.extra
        copyat_or_push!(
            k, 14,
            f(
                uprev +
                    dt *
                    (
                    a1401 * k[1] + a1406 * k[6] + a1407 * k[7] + a1408 * k[8] +
                        a1409 * k[9] + a1410 * k[10] + a1411 * k[11] + a1412 * k[12]
                ),
                p,
                t + c14 * dt
            )
        )
        copyat_or_push!(
            k, 15,
            f(
                uprev +
                    dt *
                    (
                    a1501 * k[1] + a1506 * k[6] + a1507 * k[7] + a1508 * k[8] +
                        a1509 * k[9] + a1510 * k[10] + a1511 * k[11] + a1512 * k[12] +
                        a1514 * k[14]
                ),
                p,
                t + c15 * dt
            )
        )
        copyat_or_push!(
            k, 16,
            f(
                uprev +
                    dt *
                    (
                    a1601 * k[1] + a1606 * k[6] + a1607 * k[7] + a1608 * k[8] +
                        a1609 * k[9] + a1610 * k[10] + a1611 * k[11] + a1612 * k[12] +
                        a1614 * k[14] + a1615 * k[15]
                ),
                p,
                t + c16 * dt
            )
        )
        copyat_or_push!(
            k, 17,
            f(
                uprev +
                    dt *
                    (
                    a1701 * k[1] + a1706 * k[6] + a1707 * k[7] + a1708 * k[8] +
                        a1709 * k[9] + a1710 * k[10] + a1711 * k[11] + a1712 * k[12] +
                        a1714 * k[14] + a1715 * k[15] + a1716 * k[16]
                ),
                p,
                t + c17 * dt
            )
        )
        copyat_or_push!(
            k, 18,
            f(
                uprev +
                    dt *
                    (
                    a1801 * k[1] + a1806 * k[6] + a1807 * k[7] + a1808 * k[8] +
                        a1809 * k[9] + a1810 * k[10] + a1811 * k[11] + a1812 * k[12] +
                        a1814 * k[14] + a1815 * k[15] + a1816 * k[16] + a1817 * k[17]
                ),
                p, t + c18 * dt
            )
        )
        copyat_or_push!(
            k, 19,
            f(
                uprev +
                    dt *
                    (
                    a1901 * k[1] + a1906 * k[6] + a1907 * k[7] + a1908 * k[8] +
                        a1909 * k[9] + a1910 * k[10] + a1911 * k[11] + a1912 * k[12] +
                        a1914 * k[14] + a1915 * k[15] + a1916 * k[16] + a1917 * k[17]
                ),
                p, t + c19 * dt
            )
        )
        copyat_or_push!(
            k, 20,
            f(
                uprev +
                    dt *
                    (
                    a2001 * k[1] + a2006 * k[6] + a2007 * k[7] + a2008 * k[8] +
                        a2009 * k[9] + a2010 * k[10] + a2011 * k[11] + a2012 * k[12] +
                        a2014 * k[14] + a2015 * k[15] + a2016 * k[16] + a2017 * k[17]
                ),
                p, t + c20 * dt
            )
        )
        copyat_or_push!(
            k, 21,
            f(
                uprev +
                    dt *
                    (
                    a2101 * k[1] + a2106 * k[6] + a2107 * k[7] + a2108 * k[8] +
                        a2109 * k[9] + a2110 * k[10] + a2111 * k[11] + a2112 * k[12] +
                        a2114 * k[14] + a2115 * k[15] + a2116 * k[16] + a2117 * k[17]
                ),
                p, t + c21 * dt
            )
        )
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Vern9ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    if length(k) < 10 || always_calc_begin
        @OnDemandTableauExtract Vern9Tableau T T2
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(uprev + dt * (a0201 * k[1]), p, t + c1 * dt))
        copyat_or_push!(k, 3, f(uprev + dt * (a0301 * k[1] + a0302 * k[2]), p, t + c2 * dt))
        copyat_or_push!(k, 4, f(uprev + dt * (a0401 * k[1] + a0403 * k[3]), p, t + c3 * dt))
        copyat_or_push!(
            k, 5,
            f(
                uprev + dt * (a0501 * k[1] + a0503 * k[3] + a0504 * k[4]), p,
                t + c4 * dt
            )
        )
        copyat_or_push!(
            k, 6,
            f(
                uprev + dt * (a0601 * k[1] + a0604 * k[4] + a0605 * k[5]), p,
                t + c5 * dt
            )
        )
        copyat_or_push!(
            k, 7,
            f(
                uprev +
                    dt * (a0701 * k[1] + a0704 * k[4] + a0705 * k[5] + a0706 * k[6]),
                p, t + c6 * dt
            )
        )
        copyat_or_push!(
            k, 2,
            f(
                uprev + dt * (a0801 * k[1] + a0806 * k[6] + a0807 * k[7]), p,
                t + c7 * dt
            )
        )
        copyat_or_push!(
            k, 3,
            f(
                uprev +
                    dt * (a0901 * k[1] + a0906 * k[6] + a0907 * k[7] + a0908 * k[2]),
                p, t + c8 * dt
            )
        )
        copyat_or_push!(
            k, 4,
            f(
                uprev +
                    dt *
                    (
                    a1001 * k[1] + a1006 * k[6] + a1007 * k[7] + a1008 * k[2] +
                        a1009 * k[3]
                ),
                p,
                t + c9 * dt
            )
        )
        copyat_or_push!(
            k, 5,
            f(
                uprev +
                    dt *
                    (
                    a1101 * k[1] + a1106 * k[6] + a1107 * k[7] + a1108 * k[2] +
                        a1109 * k[3] + a1110 * k[4]
                ),
                p,
                t + c10 * dt
            )
        )
        temp6 = recursivecopy(k[6])
        temp7 = recursivecopy(k[7])
        copyat_or_push!(
            k, 6,
            f(
                uprev +
                    dt *
                    (
                    a1201 * k[1] + a1206 * temp6 + a1207 * temp7 + a1208 * k[2] +
                        a1209 * k[3] + a1210 * k[4] + a1211 * k[5]
                ),
                p,
                t + c11 * dt
            )
        )
        copyat_or_push!(
            k, 7,
            f(
                uprev +
                    dt *
                    (
                    a1301 * k[1] + a1306 * temp6 + a1307 * temp7 + a1308 * k[2] +
                        a1309 * k[3] + a1310 * k[4] + a1311 * k[5] + a1312 * k[6]
                ),
                p,
                t + c12 * dt
            )
        )
        copyat_or_push!(
            k, 8,
            f(
                uprev +
                    dt *
                    (
                    a1401 * k[1] + a1406 * temp6 + a1407 * temp7 + a1408 * k[2] +
                        a1409 * k[3] + a1410 * k[4] + a1411 * k[5] + a1412 * k[6] +
                        a1413 * k[7]
                ),
                p,
                t + c13 * dt
            )
        )
        copyat_or_push!(
            k, 9,
            f(
                uprev +
                    dt *
                    (
                    a1501 * k[1] + a1506 * temp6 + a1507 * temp7 + a1508 * k[2] +
                        a1509 * k[3] + a1510 * k[4] + a1511 * k[5] + a1512 * k[6] +
                        a1513 * k[7] + a1514 * k[8]
                ),
                p,
                t + dt
            )
        )
        copyat_or_push!(
            k, 10,
            f(
                uprev +
                    dt *
                    (
                    a1601 * k[1] + a1606 * temp6 + a1607 * temp7 + a1608 * k[2] +
                        a1609 * k[3] + a1610 * k[4] + a1611 * k[5] + a1612 * k[6] +
                        a1613 * k[7]
                ),
                p,
                t + dt
            )
        )
    end
    if (allow_calc_end && length(k) < 20) || force_calc_end # Have not added the extra stages yet
        @OnDemandTableauExtract Vern9ExtraStages T T2
        copyat_or_push!(
            k, 11,
            f(
                uprev +
                    dt *
                    (
                    a1701 * k[1] + a1708 * k[2] + a1709 * k[3] + a1710 * k[4] +
                        a1711 * k[5] + a1712 * k[6] + a1713 * k[7] + a1714 * k[8] +
                        a1715 * k[9]
                ),
                p,
                t + c17 * dt
            )
        )
        copyat_or_push!(
            k, 12,
            f(
                uprev +
                    dt *
                    (
                    a1801 * k[1] + a1808 * k[2] + a1809 * k[3] + a1810 * k[4] +
                        a1811 * k[5] + a1812 * k[6] + a1813 * k[7] + a1814 * k[8] +
                        a1815 * k[9] + a1817 * k[11]
                ),
                p,
                t + c18 * dt
            )
        )
        copyat_or_push!(
            k, 13,
            f(
                uprev +
                    dt *
                    (
                    a1901 * k[1] + a1908 * k[2] + a1909 * k[3] + a1910 * k[4] +
                        a1911 * k[5] + a1912 * k[6] + a1913 * k[7] + a1914 * k[8] +
                        a1915 * k[9] + a1917 * k[11] + a1918 * k[12]
                ),
                p,
                t + c19 * dt
            )
        )
        copyat_or_push!(
            k, 14,
            f(
                uprev +
                    dt *
                    (
                    a2001 * k[1] + a2008 * k[2] + a2009 * k[3] + a2010 * k[4] +
                        a2011 * k[5] + a2012 * k[6] + a2013 * k[7] + a2014 * k[8] +
                        a2015 * k[9] + a2017 * k[11] + a2018 * k[12] + a2019 * k[13]
                ),
                p,
                t + c20 * dt
            )
        )
        copyat_or_push!(
            k, 15,
            f(
                uprev +
                    dt *
                    (
                    a2101 * k[1] + a2108 * k[2] + a2109 * k[3] + a2110 * k[4] +
                        a2111 * k[5] + a2112 * k[6] + a2113 * k[7] + a2114 * k[8] +
                        a2115 * k[9] + a2117 * k[11] + a2118 * k[12] + a2119 * k[13] +
                        a2120 * k[14]
                ),
                p,
                t + c21 * dt
            )
        )
        copyat_or_push!(
            k, 16,
            f(
                uprev +
                    dt *
                    (
                    a2201 * k[1] + a2208 * k[2] + a2209 * k[3] + a2210 * k[4] +
                        a2211 * k[5] + a2212 * k[6] + a2213 * k[7] + a2214 * k[8] +
                        a2215 * k[9] + a2217 * k[11] + a2218 * k[12] + a2219 * k[13] +
                        a2220 * k[14] + a2221 * k[15]
                ),
                p,
                t + c22 * dt
            )
        )
        copyat_or_push!(
            k, 17,
            f(
                uprev +
                    dt *
                    (
                    a2301 * k[1] + a2308 * k[2] + a2309 * k[3] + a2310 * k[4] +
                        a2311 * k[5] + a2312 * k[6] + a2313 * k[7] + a2314 * k[8] +
                        a2315 * k[9] + a2317 * k[11] + a2318 * k[12] + a2319 * k[13] +
                        a2320 * k[14] + a2321 * k[15]
                ),
                p,
                t + c23 * dt
            )
        )
        copyat_or_push!(
            k, 18,
            f(
                uprev +
                    dt *
                    (
                    a2401 * k[1] + a2408 * k[2] + a2409 * k[3] + a2410 * k[4] +
                        a2411 * k[5] + a2412 * k[6] + a2413 * k[7] + a2414 * k[8] +
                        a2415 * k[9] + a2417 * k[11] + a2418 * k[12] + a2419 * k[13] +
                        a2420 * k[14] + a2421 * k[15]
                ),
                p,
                t + c24 * dt
            )
        )
        copyat_or_push!(
            k, 19,
            f(
                uprev +
                    dt *
                    (
                    a2501 * k[1] + a2508 * k[2] + a2509 * k[3] + a2510 * k[4] +
                        a2511 * k[5] + a2512 * k[6] + a2513 * k[7] + a2514 * k[8] +
                        a2515 * k[9] + a2517 * k[11] + a2518 * k[12] + a2519 * k[13] +
                        a2520 * k[14] + a2521 * k[15]
                ),
                p,
                t + c25 * dt
            )
        )
        copyat_or_push!(
            k, 20,
            f(
                uprev +
                    dt *
                    (
                    a2601 * k[1] + a2608 * k[2] + a2609 * k[3] + a2610 * k[4] +
                        a2611 * k[5] + a2612 * k[6] + a2613 * k[7] + a2614 * k[8] +
                        a2615 * k[9] + a2617 * k[11] + a2618 * k[12] + a2619 * k[13] +
                        a2620 * k[14] + a2621 * k[15]
                ),
                p,
                t + c26 * dt
            )
        )
    end
    nothing
end
