.SECONDEXPANSION:

define directories

$1: | $$$$(dir $$$$@)
	mkdir $$@

endef

#${directories}: | $$(dir $$@)
#	mkdir $@

./:
	@# Nothing to be done.

# vim: ft=make
