
ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: build/debug/atf build/release/atf

build/debug/atf: build/debug atf.tpl 
	@admb -g atf.tpl 
	@$(COPY) atf build/debug/
	@$(RM) atf

build/release/atf: build/release atf.tpl 
	@admb -g atf.tpl 
	@$(COPY) atf build/release/
	@$(RM) atf

build/debug:
	@echo creating debug directory
	@mkdir -p build/debug

build/release:
	@echo creating release directory
	@mkdir -p build/release

clean:
	@$(RM) build
	@$(RM) atf.cpp
	@$(RM) atf.htp
	@$(RM) atf.obj

