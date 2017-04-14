ifeq ($(KAM),TRUE)
  EXEC = kam
else
  EXEC = atf
endif
ARGS = -nox -iprint 50 

ifeq ($(SAFE),FALSE)
	DIST = ../src/build/release/
else
	DIST = ../src/build/debug/
endif

ifdef ComSpec
    RM=del /F /Q
else
    RM=rm -rf
endif

all: run

$(EXEC): $(DIST)$(EXEC).tpl
	ln -sf $(DIST)$(EXEC) $@

$(DIST)$(EXEC).tpl:
	$(MAKE) --directory=../src

run: $(EXEC)
	./$(EXEC) $(ARGS)

plot:
	R CMD BATCH ..\R\ATF_plotting.R

both: $(EXEC)
	./$(EXEC) $(ARGS)
	R CMD BATCH plot.R

sim: $(EXEC)
	./$(EXEC) -ainp atf.par -sim 123

clean:
	@$(RM) -f $(EXEC)
	@$(RM) -f admodel.*
	@$(RM) -f $(EXEC)  $(EXEC).[brces]* $(EXEC).*[0123456789] *.rpt *.log variance

clobber: clean
	@$(RM) -f *.par *_in.* *.rep

