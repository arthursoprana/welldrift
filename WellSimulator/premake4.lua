project "WellSimulator"

kind "ConsoleApp"
targetdir "./Bin"
	
configuration {}

includedirs {
	"./",
	"./Solver",
	"./WellSim/include",
	"./WinApi",
	"../libgtkgraph",
}
libdirs {
}
links {
	"gtkgraph",
	"glib-2.0",
	"gtk-win32-2.0",
	"gdk-win32-2.0",
	"gobject-2.0",
	"gdk_pixbuf-2.0",
	"gthread-2.0",
	"gmodule-2.0",
	"pango-1.0",
	"atk-1.0",
	"zdll",	
}

files {
	"**.cpp",
	"**.c",
	"**.cc",
	"**.h",
	"**.hpp"
}
