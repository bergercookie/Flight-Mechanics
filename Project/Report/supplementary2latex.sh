# Make nomenclature
makeindex report.nlo -s nomencl.ist -o report.nls -t report.nlg
makeindex report.nlo -s nomencl.ist -o report.nls -t report.nlg

# Render the bibliography section
bibtex report
bibtex report

# Run xelatex compuler
xelatex report

# Once more for the bibliography, contents to appear correctly
echo "******************************************************"
echo "*                    ONCE AGAIN                      *"
echo "******************************************************"
xelatex report

