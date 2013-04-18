tests = test1

.PHONY += selfcheck $(tests)

selfcheck : $(tests)
	@echo Done testing!

test1 : include_modelHeader:=$(selfcheck_dir)test1/
test% : model
	$(SHELL) $(selfcheck_dir)runselfcheck.sh $@
