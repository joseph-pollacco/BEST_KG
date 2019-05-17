
using Pkg
function PACKAGES(;Option_PackageUpdate = false)


    function PACKAGE_MANAGER(Package)
		"""ONLY ADD PACKAGE IF NOT AVAILABLE"""
			if  haskey(Pkg.installed(), Package) == false
				println("Adding $Package package because not available...")
				Pkg.add(Package)
				println("$Package package added")
			end
		end # PACKAGE_MANAGER

#    setprotocol! # git config --global url."https://github.com/".insteadOf git://github.com/
#     ENV["PYTHON"]="" # Setting python environment
    PACKAGE_MANAGER("NLopt")
    # PACKAGE_MANAGER("PyCall")
    PACKAGE_MANAGER("CSV")
    PACKAGE_MANAGER("SpecialFunctions")
    PACKAGE_MANAGER("QuadGK")
    # PACKAGE_MANAGER("Plots")
    # PACKAGE_MANAGER("TikzGraphs")
    # Pkg.checkout("Plots")
    # PACKAGE_MANAGER("PyPlot")
    PACKAGE_MANAGER("PGFPlots")
    # PACKAGE_MANAGER("GR")
    # PACKAGE_MANAGER("UnicodePlots")
    # PACKAGE_MANAGER("StatPlots")
    # PACKAGE_MANAGER("PlotlyJS")
    # PACKAGE_MANAGER("PlotRecipes")
    PACKAGE_MANAGER("BlackBoxOptim")
    PACKAGE_MANAGER("ExcelReaders")
    # PACKAGE_MANAGER("Cubature")
    PACKAGE_MANAGER("Optim")
    # PACKAGE_MANAGER("SymEngine")
    # PACKAGE_MANAGER("Roots")
    # PACKAGE_MANAGER("QuantEcon")

    if Option_PackageUpdate
        println("Updating metdata...")
        Pkg.update()
    end

end

 PACKAGES(Option_PackageUpdate = true)
