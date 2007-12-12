.onAttach <- function(libname, pkgname)
#	Add HowToPLW to Windows menu
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) winMenuAddItem("Vignettes","plw","HowToPLW()")
}
