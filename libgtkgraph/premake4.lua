project "gtkgraph"

kind "StaticLib"
targetdir "../lib"

	
configuration {}
includedirs {
		"./",
		"../gtk+/include/",
		"../gtk+/include/atk-1.0",
		"../gtk+/include/glib-2.0",
		"../gtk+/include/gtk-2.0",
		"../gtk+/include/gdk-pixbuf-2.0",
		"../gtk+/include/cairo",
		"../gtk+/include/pango-1.0",
		"../gtk+/lib/gtk-2.0/include",
		"../gtk+/lib/glib-2.0/include",
	}
	libdirs {
		"../gtk+/lib",	
	}
	links{
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
	"**.cc",
	"**.c",
	"**.h"
}