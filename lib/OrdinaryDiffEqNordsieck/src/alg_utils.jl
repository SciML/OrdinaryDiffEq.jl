alg_order(alg::AN5) = 5
alg_order(alg::JVODE) = 1  #dummy value

alg_adaptive_order(alg::AN5) = 5

qsteady_max_default(alg::AN5) = 3 // 2
qsteady_max_default(alg::JVODE) = 3 // 2

get_current_alg_order(alg::JVODE, cache) = get_current_adaptive_order(alg, cache)

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function legacy_default_controller(alg::JVODE, args...)
            return DummyController()
        end

        function default_controller_v7(QT, alg::JVODE, args...)
            return DummyController()
        end
    end
else
    @eval begin
        function OrdinaryDiffEqCore.default_controller(alg::JVODE, args...)
            return DummyController()
        end
    end
end
