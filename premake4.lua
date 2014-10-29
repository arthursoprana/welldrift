solution "WellDrift"

configurations {"Release", "Debug"}
configuration "Release"
	flags {"StaticRuntime", "NoMinimalRebuild", "NoFramePointer", "OptimizeSpeed"}
configuration "Debug"
	flags {"StaticRuntime", "Symbols","NoMinimalRebuild", "NoEditAndContinue"}

platforms {"x32","x64"}

configuration {"x64", "Debug"}
		targetsuffix "x64_d"
		
configuration {"x64", "Release"}
		targetsuffix "x64"		

configuration {"Windows"}
	defines {"WIN32", "_CONSOLE", "_CRT_SECURE_NO_WARNINGS"}
	
configuration {"Debug"}
	defines {"_DEBUG"}
configuration {"Release"}
	defines {"NDEBUG"}
configuration{}

language "C++"

location ("./Solution")


includedirs {
	"./",	
	"../boost_1_56_0",
	"../gtk+/include",
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
	"./lib",
	"../gtk+/lib",
}

include "./WellSimulator"
include "./libgtkgraph"







