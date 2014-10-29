solution "WellDrift"

configurations {"Release", "Debug"}
configuration "Release"
	flags {"StaticRuntime", "NoMinimalRebuild", "Unicode", "NoFramePointer", "OptimizeSpeed"}
configuration "Debug"
	flags {"StaticRuntime", "Symbols","NoMinimalRebuild", "Unicode", "NoEditAndContinue"}

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
	"../../boost_1_56_0",
}
libdirs {
	"./lib",	
}

include "./WellSimulator"
include "./libgtkgraph"







