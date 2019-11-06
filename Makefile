DOCS=index research teaching resume

HDOCS=$(addsuffix .html, $(DOCS))

.PHONY : docs
docs : $(HDOCS)

.PHONY : update
update : $(PHDOCS)
	@echo -n 'Copying to server...'
	git add *
	git commit -m "website update"
	git push origin master
	@echo ' done.'

%.html : %.jemdoc MENU
	@echo -n 'Running jemdoc.'
	python2.7 jemdoc.py -o $@ $<
	python2.7 jemdoc.py -c mysite.conf -o index.html index.jemdoc

.PHONY : clean
clean :
	-rm -f *.html
