import requests
from bs4 import BeautifulSoup


# params
bed_file = '/home/alabarga/NavarraBiomed/analysis/parkinson/results/selected.cpgsPDvsCTRL.bed'


# GREAT
url = 'http://great.stanford.edu/public/cgi-bin/greatWeb.php'
headers = {'user-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/64.0.3282.119 Chrome/64.0.3282.119 Safari/537.36'}
params = { 'species': 'hg19', 
          'rule': 'basalPlusExt',
          'span': '1000.0',
          'upstream': '5.0',
          'downstream': '1.0',
          'twoDistance': '1000.0',
          'oneDistance': '1000.0',
          'includeCuratedRegDoms': '1',
          'bgChoice': 'wholeGenome',
          'fgChoice': 'file',
        }

files = {'fgFile': open(bed_file, 'r')}

# launch job
r = requests.post(url, headers=headers, data=params, files=files)

if r.status_code:
    soup = BeautifulSoup(r.text, 'lxml')
    jobId = soup.find('div', attrs={'class': 'job_desc_info'}).text



    # GREAT
    downloadUrl = 'http://great.stanford.edu/public/cgi-bin/downloadAllTSV.php'
    headers = {'user-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/64.0.3282.119 Chrome/64.0.3282.119 Safari/537.36'}
    params = { 'outputDir': '/scratch/great/tmp/results/%s/' % (jobId), 
              'outputPath': '/tmp/results/%s/' % (jobId),
              'genomeAssemblyDir': '/webapp/great/public-3.0.0-ontologies/hg19/20140326',
              'species': 'hg19',          
              'ontoName': '',
              'ontoList': 'GOMolecularFunction@1-Inf,GOBiologicalProcess@1-Inf,GOCellularComponent@1-Inf,MGIPhenotype@1-Inf,HumanPhenotypeOntology@1-Inf,OsborneAnnotatedDiseaseOntology@1-Inf,MSigDBGeneSetsCancerNeighborhood@1-Inf,PlacentaDisorders@1-Inf,panther@1-Inf,BioCycPathway@1-Inf,MSigDBGeneSetsCanonicalPathway@1-Inf,MGIExpressionDetected@1-Inf,MSigDBGeneSetsPerturbation@1-Inf,MSigDBGeneSetsPromoterMotifs@1-Inf,MSigDBGeneSetsMicroRNAMotifs@1-Inf,interpro@1-Inf,treefam@1-Inf,geneFamilies@1-Inf,EnsemblGenes@1-Inf,MSigDBGeneSetsOncogenicSignatures@1-Inf,MSigDBGeneSetsImmunologicSignatures@1-Inf',
              'binom': 'true',
            }

    r2 = requests.post(downloadUrl, headers=headers, data=params)

    f = open("results.tsv",'w')
    f.write(r2.text)
    f.close()
