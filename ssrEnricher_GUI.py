#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Date: 2020.05.09
# Author: Wu Qing
# SSREnricher_GUI
# Version: V1.1
# This software can automatically enrich polymorphic SSRs with transcripts.
# -----------------------------------------------------------------------------
from Bio import SeqIO
import os
import re
from collections import defaultdict
import sys
from tkinter import *
from tkinter.ttk import *
from tkinter.messagebox import *
import tkinter.filedialog
import time
import subprocess

def canRunCommand(command, qAllowStderr = False, qPrint = True):
    if qPrint: sys.stdout.write("Test can run \"%s\"" % command)
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=os.environ.copy())
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
        if qPrint: print(" - ok")
        return True
    else:
        if qPrint: print(" - failed")
        print("\nstdout:")        
        for l in stdout: print(l)
        print("\nstderr:")        
        for l in stderr: print(l)
        return False

def checkPerl():
    if canRunCommand('perl -h'):
        return True
    else:
        print("ERROR: Cannot run perl")
        print("Please check perl is installed and that the executables are in the system path\n")
        showinfo(
        title='info',
        message="ERROR: Cannot run perl\n\nPlease check perl is installed and that the executables are in the system path\n")
        raise
        return False

def checkCdhit():
    if canRunCommand('cd-hit-est -h'):
        return True
    else:
        print("ERROR: Cannot run cd-hit-est")
        print("Please check cd-hit is installed and that the executables are in the system path\n")
        showinfo(
        title='info',
        message="ERROR: Cannot run cd-hit-est\n\nPlease check cd-hit is installed and that the executables are in the system path\n")
        raise
        return False

def checkMuscle():
    if canRunCommand('muscle -version'):
        return True
    else:
        print("ERROR: Cannot run muscle")
        print("Please check muscle is installed and that the executables are in the system path\n")
        showinfo(
        title='info',
        message="ERROR: Cannot run muscle\n\nPlease check muscle is installed and that the executables are in the system path\n")
        raise
        return False

def getSamples(n):
    namelist = []
    for fn in os.listdir(os.getcwd()):
        if os.path.splitext(fn)[1] == n:
            fn = ''.join(os.path.splitext(fn))
            namelist.append(fn)
    return namelist


def fixGeneId():
    os.system('rm -rf WorkingDrictory && mkdir WorkingDrictory')
    os.chdir('./WorkingDrictory')
    i = 1
    for sample in samplesFile:
        fas = SeqIO.parse(open(sample), 'fasta')
        for fa in fas:
            a = os.getcwd()
            open('T' + str(i) + '.fa', 'a').write('>T%s.%s\n%s\n' %(i, '_'.join(fa.id.split('_')[1:]), str(fa.seq)))
        i += 1


def callMisa(samples):
    for i in samples:
        os.system('perl ./misa.pl ' + i)


def getSsrSeq(samples):
    # get SSR sequcens from fasta by misa result ; combine all group
    #samples = ['T1.fa', 'T3.fa', 'T2.fa']

    for sample in samples:
        ids = []
        faD = SeqIO.to_dict(SeqIO.parse(open(sample), 'fasta'))
        for la in open(sample + '.misa'):
            if 'ID' not in la:
                aL = la.split('\t')
                ln = len(faD[aL[0]].seq)
                if aL[2] != 'p1' and aL[5] != 1 and aL[6] != ln:
                    ids.append(aL[0])
        comp_ids = []
        for lb in open(sample + '.misa'):
            bL = lb.strip().split('\t')
            if 'c' in bL[2]:
                comp_ids.append(bL[0])
        fas = SeqIO.parse(open(sample), 'fasta')
        for fa in fas:
            if fa.id in ids and fa.id not in comp_ids:
                open(
                    'all.ssr.fa', 'a').write(
                    '>%s\n%s\n' %
                    (fa.id, str(
                        fa.seq)))


def callCdHit():
    # check cd-hit depende and search the cluster sequences by CD-HIT-EST
    os.system('cd-hit-est -i all.ssr.fa -o cd-hit')


def reformSSR(samples):
    #samples = ['T1.fa', 'T2.fa', 'T3.fa']
    faD = SeqIO.to_dict(SeqIO.parse(open('all.ssr.fa'), 'fasta'))
    baseD = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for sample in samples:
        for line in open(sample + '.misa'):
            lst = line.strip().split('\t')
            if lst[0] in faD:
                if 'c' not in lst[2] and lst[2] != 'p1' and 'ID' not in line and int(lst[-2]) > 50 and int(lst[-1]) < len(
                    faD[lst[0]].seq) - 50:
                    ma = re.findall(r'\((.+)\)', lst[3])[0]
                    maL = list(ma)
                    rc = ''
                    for base in maL:
                        if base in 'ACGT':
                            rc += baseD[base]
                    rc = rc[::-1]
                    ss = sorted([ma, rc])
                    ss = '(' + '/'.join(ss) + ')'
                    newline = re.sub(r'\(([ACGT]+)\)', ss, line)
                    open(sample + '.ssr.reformed', 'a').write(str(newline))


def getReverseSeq():
    sD = {}
    for la in open('cd-hit.clstr'):
        if 'at' in la:
            id = re.findall(r'>(.+)\.\.\.', la)[0]
            strand = re.findall(r'([+-])+\/', la)[0]
            sD[id] = strand
        if '*' in la:
            id = re.findall(r'>(.+)\.\.\.', la)[0]
            sD[id] = '+'
    fas = SeqIO.parse(open('all.ssr.fa'), 'fasta')
    for fa in fas:
        if fa.id in sD:
            if sD[fa.id] == '+':
                open('plus.ssr.fa', 'a').write(
                    '>' + str(fa.id) + '\n' + str(fa.seq) + '\n')
            if sD[fa.id] == '-':
                seq = fa.seq.reverse_complement()
                open('plus.ssr.fa', 'a').write(
                    '>' + str(fa.id) + '\n' + str(seq) + '\n')

def enrichSSR():
    def getD(ssr):
        s = []
        for la in open(ssr):
            if 'ID' not in la:
                aL = la.strip().split('\t')
                ma = re.findall(r'\(.+\)\d+', aL[3])
                s.append((aL[0], ma[0]))
        d = defaultdict(set)
        for k, v in s:
            d[k].add(v)
        return d

    allD = {}

    for i in getSamples('.reformed'):
        t = getD(i)
        allD.update(t)

    page = open('cd-hit.clstr').read()
    clusters = re.findall('(.+?)>Cluster', page, re.S)
    for cluster in clusters:
        trans = re.findall(r'T\d.+c\d+_g\d+_i\d+', cluster)
        if len(trans) > 1:
            tt = []
            ss = []
            for tran in trans:
                if tran in allD:
                    tt.append(str(tran) + ':' + str(list(allD[tran])))
                    ss += list(allD[tran])
            if len(tt) > 1:
                ma = re.findall(r'\)\d+\'', str(ss))
                ma = set(ma)
                mas = re.findall(r'\((.+?)\)', str(set(ss)))
                ssr = ''
                for mm in list(set(mas)):
                    if mas.count(mm) > 1:
                        ssr = mm

                if len(ma) > 1 and len(mas) > len(set(mas)):
                    ttt = []
                    for t in tt:
                        if ssr in t:
                            ttt.append(t)
                    na = str(ttt).lstrip('[').rstrip(']').strip(
                        '"').replace('", "', '\t')
                    open('enrichment.SSRs', 'a').write(str(na) + '\n')

    for ll in open('enrichment.SSRs'):
        ma = re.findall(r'T\d\.[a-zA-Z0-9]+_c\d+_g\d+_i\d+', ll)
        ma = set(ma)
        if len(ma) > 1:
            open('enrichment.SSRs.txt', 'a').write(str(ll))
    
    #global resultDir
    resultDir = time.strftime("Result_%Y-%m-%d_%H-%M", time.localtime())
    os.system('mkdir ../' + resultDir)

    with open('enrichment.SSRs.txt') as f:
        for i in f:
            a = re.findall(r'\[(.+?)\]',i)
            if len(set(a)) !=1:
                sdic = {}
                for k in a:
                    q = re.findall(r'\((.+?)\)', k)
                    q = '-'.join(q)
                    v = re.findall(r'\d+',k)
                    v = '-'.join(v)
                    l =[]
                    if q not in sdic:
                        l.append(v)
                        sdic[q] = l
                    else:
                        sdic[q].append(v)
                if len(sdic) == 1 :
                    c =  len(set(sdic[q]))
                    if c != 1:
                        open('SSR_Result.txt','a').write(i)
                os.system('cp SSR_Result.txt ../' + resultDir)


    faD = SeqIO.to_dict(SeqIO.parse(open('plus.ssr.fa'), 'fasta'))
    n = 0
    for la in open('SSR_Result.txt'):
        aL = la.strip().split('\t')
        for it in aL:
            id = it.split(':')[0]
            open('seq' + str(n), 'a').write('>' + \
                 str(id) + '\n' + str(faD[id].seq) + '\n')
        os.system(
            'muscle -msf -in seq' +
            str(n) +
            ' -out seq' +
            str(n) +
            '.muscle')
        n += 1
    os.system('mkdir ../'+ resultDir + '/muscle_Result ; mv seq* ../'+ resultDir + '/muscle_Result; rm enrichment.SSRs*')
    return resultDir

def getFinal(resultDir):
    faDic = SeqIO.to_dict(SeqIO.parse('all.ssr.fa', 'fasta'))
    ssrDic = {}
    os.system('cat *.reformed > all.reformed')
    with open ('./all.reformed', 'r') as f:
        for line in f:
            a = line.split('\t')[0]
            b = line.split('\t')[5:]
            if a not in ssrDic:
                ssrDic[a] = '-'.join(b)[:-1]
            else:
                ssrDic[a] += '\t' + '-'.join(b)[:-1]

                

    with open ('./SSR_Result.txt', 'r') as f:
        for line in f:
            with open('Fasta_Result.txt', 'a') as f:
                f.write('\n' +'-'*80 + '\n' + line + '\n')

            a = line.split('\t')
            for i in a :
                id = i.split(':')[0]
                ssr = i.split(':')[1]

                with open('Fasta_Result.txt', 'a') as f:
                    f.write(i + '\t' + ssrDic[id] + '\n' + str(faDic[id].seq + '\n\n'))
    
    os.system('cp Fasta_Result.txt ../' + resultDir)

def selectCom_Cmd():
    global samplesFile
    samplesFile = tkinter.filedialog.askopenfilenames()
    fileNa = ''
    for i in samplesFile:
        fileNa = fileNa + i.split("/")[-1] + '\n'
    SelectLabel.config(text="You choosed：\n" + fileNa)

def creatMisa():
    with open('./misa.pl','w') as f:
        f.write('''#!/usr/bin/perl -w\n# Author: Thomas Thiel\n# Program name: misa.pl\n\n###_______________________________________________________________________________\n###\n###Program name: misa.pl\n###Author:       Thomas Thiel\n###Release date: 14/12/01 (version 1.0)\n###\n###_______________________________________________________________________________\n###\n## _______________________________________________________________________________\n##\n## DESCRIPTION: Tool for the identification and localization of\n##              (I)  perfect microsatellites as well as\n##              (II) compound microsatellites (two individual microsatellites,\n##                   disrupted by a certain number of bases)\n##\n## SYNTAX:   misa.pl <FASTA file>\n##\n##    <FASTAfile>    Single file in FASTA format containing the sequence(s).\n##\n##    In order to specify the search criteria, an additional file containing\n##    the microsatellite search parameters is required named "misa.ini", which\n##    has the following structure:\n##      (a) Following a text string beginning with \'def\', pairs of numbers are\n##          expected, whereas the first number defines the unit size and the\n##          second number the lower threshold of repeats for that specific unit.\n##      (b) Following a text string beginning with \'int\' a single number defines\n##          the maximal number of bases between two adjacent microsatellites in\n##          order to specify the compound microsatellite type.\n##    Example:\n##      definition(unit_size,min_repeats):          1-10 2-6 3-5 4-5 5-5 6-5\n##      interruptions(max_difference_for_2_SSRs):   100\n##\n## EXAMPLE: misa.pl seqs.fasta\n##\n## _______________________________________________________________________________\n##\n\n\n#§§§§§ DECLARATION §§§§§#\n\n# Check for arguments. If none display syntax #\n\nif (@ARGV == 0)\n  {\n  open (IN,"<$0");\n  while (<IN>) {if (/^\\#\\# (.*)/) {$message .= "$1\\n"}};\n  close (IN);\n  die $message;\n  };\n\n# Check if help is required #\n\nif ($ARGV[0] =~ /-help/i)\n  {\n  open (IN,"<$0");\n  while (<IN>) {if (/^\\#\\#\\#(.*)/) {$message .= "$1\\n"}};\n  close (IN);\n  die $message;\n  };\n\n# Open FASTA file #\n\nopen (IN,"<$ARGV[0]") || die ("\\nError: FASTA file doesn\'t exist !\\n\\n");\nopen (OUT,">$ARGV[0].misa");\nprint OUT "ID\\tSSR nr.\\tSSR type\\tSSR\\tsize\\tstart\\tend\\n";\n\n# Reading arguments #\n\nopen (SPECS,"misa.ini") || die ("\\nError: Specifications file doesn\'t exist !\\n\\n");\nmy %typrep;\nmy $amb = 0;\nwhile (<SPECS>)\n   {\n   %typrep = $1 =~ /(\\d+)/gi if (/^def\\S*\\s+(.*)/i);\n   if (/^int\\S*\\s+(\\d+)/i) {$amb = $1}\n   };\nmy @typ = sort { $a <=> $b } keys %typrep;\n\n\n#§§§§§ CORE §§§§§#\n\n$/ = ">";\nmy $max_repeats = 1; #count repeats\nmy $min_repeats = 1000; #count repeats\nmy (%count_motif,%count_class); #count\nmy ($number_sequences,$size_sequences,%ssr_containing_seqs); #stores number and size of all sequences examined\nmy $ssr_in_compound = 0;\nmy ($id,$seq);\nwhile (<IN>)\n  {\n  next unless (($id,$seq) = /(.*?)\\n(.*)/s);\n  my ($nr,%start,@order,%end,%motif,%repeats); # store info of all SSRs from each sequence\n  $seq =~ s/[\\d\\s>]//g; #remove digits, spaces, line breaks,...\n  $id =~ s/^\\s*//g; $id =~ s/\\s*$//g;$id =~ s/\\s/_/g; #replace whitespace with "_"\n  $number_sequences++;\n  $size_sequences += length $seq;\n  for ($i=0; $i < scalar(@typ); $i++) #check each motif class\n    {\n    my $motiflen = $typ[$i];\n    my $minreps = $typrep{$typ[$i]} - 1;\n    if ($min_repeats > $typrep{$typ[$i]}) {$min_repeats = $typrep{$typ[$i]}}; #count repeats\n    my $search = "(([acgt]{$motiflen})\\\\2{$minreps,})";\n    while ( $seq =~ /$search/ig ) #scan whole sequence for that class\n      {\n      my $motif = uc $2;\n      my $redundant; #reject false type motifs [e.g. (TT)6 or (ACAC)5]\n      for ($j = $motiflen - 1; $j > 0; $j--)\n        {\n        my $redmotif = "([ACGT]{$j})\\\\1{".($motiflen/$j-1)."}";\n        $redundant = 1 if ( $motif =~ /$redmotif/ )\n        };\n      next if $redundant;\n      $motif{++$nr} = $motif;\n      my $ssr = uc $1;\n      $repeats{$nr} = length($ssr) / $motiflen;\n      $end{$nr} = pos($seq);\n      $start{$nr} = $end{$nr} - length($ssr) + 1;\n      # count repeats\n      $count_motifs{$motif{$nr}}++; #counts occurrence of individual motifs\n      $motif{$nr}->{$repeats{$nr}}++; #counts occurrence of specific SSR in its appearing repeat\n      $count_class{$typ[$i]}++; #counts occurrence in each motif class\n      if ($max_repeats < $repeats{$nr}) {$max_repeats = $repeats{$nr}};\n      };\n    };\n  next if (!$nr); #no SSRs\n  $ssr_containing_seqs{$nr}++;\n  @order = sort { $start{$a} <=> $start{$b} } keys %start; #put SSRs in right order\n  $i = 0;\n  my $count_seq; #counts\n  my ($start,$end,$ssrseq,$ssrtype,$size);\n  while ($i < $nr)\n    {\n    my $space = $amb + 1;\n    if (!$order[$i+1]) #last or only SSR\n      {\n      $count_seq++;\n      my $motiflen = length ($motif{$order[$i]});\n      $ssrtype = "p".$motiflen;\n      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";\n      $start = $start{$order[$i]}; $end = $end{$order[$i++]};\n      next\n      };\n    if (($start{$order[$i+1]} - $end{$order[$i]}) > $space)\n      {\n      $count_seq++;\n      my $motiflen = length ($motif{$order[$i]});\n      $ssrtype = "p".$motiflen;\n      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";\n      $start = $start{$order[$i]}; $end = $end{$order[$i++]};\n      next\n      };\n    my ($interssr);\n    if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)\n      {\n      $count_seq++; $ssr_in_compound++;\n      $ssrtype = \'c*\';\n      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}($motif{$order[$i+1]})$repeats{$order[$i+1]}*";\n      $start = $start{$order[$i]}; $end = $end{$order[$i+1]}\n      }\n    else\n      {\n      $count_seq++; $ssr_in_compound++;\n      $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);\n      $ssrtype = \'c\';\n      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";\n      $start = $start{$order[$i]};  $end = $end{$order[$i+1]};\n      #$space -= length $interssr\n      };\n    while ($order[++$i + 1] and (($start{$order[$i+1]} - $end{$order[$i]}) <= $space))\n      {\n      if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)\n        {\n        $ssr_in_compound++;\n        $ssrseq .= "($motif{$order[$i+1]})$repeats{$order[$i+1]}*";\n        $ssrtype = \'c*\';\n        $end = $end{$order[$i+1]}\n        }\n      else\n        {\n        $ssr_in_compound++;\n        $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);\n        $ssrseq .= "$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";\n        $end = $end{$order[$i+1]};\n        #$space -= length $interssr\n        }\n      };\n    $i++;\n    }\n  continue\n    {\n    print OUT "$id\\t$count_seq\\t$ssrtype\\t$ssrseq\\t",($end - $start + 1),"\\t$start\\t$end\\n"\n    };\n  };\n\nclose (OUT);\nopen (OUT,">$ARGV[0].statistics");\n\n#§§§§§ INFO §§§§§#\n\n#§§§ Specifications §§§#\nprint OUT "Specifications\\n==============\\n\\nSequence source file: \\"$ARGV[0]\\"\\n\\nDefinement of microsatellites (unit size / minimum number of repeats):\\n";\nfor ($i = 0; $i < scalar (@typ); $i++) {print OUT "($typ[$i]/$typrep{$typ[$i]}) "};print OUT "\\n";\nif ($amb > 0) {print OUT "\\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  $amb\\n"};\nprint OUT "\\n\\n\\n";\n\n#§§§ OCCURRENCE OF SSRs §§§#\n\n#small calculations\nmy @ssr_containing_seqs = values %ssr_containing_seqs;\nmy $ssr_containing_seqs = 0;\nfor ($i = 0; $i < scalar (@ssr_containing_seqs); $i++) {$ssr_containing_seqs += $ssr_containing_seqs[$i]};\nmy @count_motifs = sort {length ($a) <=> length ($b) || $a cmp $b } keys %count_motifs;\nmy @count_class = sort { $a <=> $b } keys %count_class;\nfor ($i = 0; $i < scalar (@count_class); $i++) {$total += $count_class{$count_class[$i]}};\n\n#§§§ Overview §§§#\nprint OUT "RESULTS OF MICROSATELLITE SEARCH\\n================================\\n\\n";\nprint OUT "Total number of sequences examined:              $number_sequences\\n";\nprint OUT "Total size of examined sequences (bp):           $size_sequences\\n";\nprint OUT "Total number of identified SSRs:                 $total\\n";\nprint OUT "Number of SSR containing sequences:              $ssr_containing_seqs\\n";\nprint OUT "Number of sequences containing more than 1 SSR:  ",$ssr_containing_seqs - ($ssr_containing_seqs{1} || 0),"\\n";\nprint OUT "Number of SSRs present in compound formation:    $ssr_in_compound\\n\\n\\n";\n\n#§§§ Frequency of SSR classes §§§#\nprint OUT "Distribution to different repeat type classes\\n---------------------------------------------\\n\\n";\nprint OUT "Unit size\\tNumber of SSRs\\n";\nmy $total = undef;\nfor ($i = 0; $i < scalar (@count_class); $i++) {print OUT "$count_class[$i]\\t$count_class{$count_class[$i]}\\n"};\nprint OUT "\\n";\n\n#§§§ Frequency of SSRs: per motif and number of repeats §§§#\nprint OUT "Frequency of identified SSR motifs\\n----------------------------------\\n\\nRepeats";\nfor ($i = $min_repeats;$i <= $max_repeats; $i++) {print OUT "\\t$i"};\nprint OUT "\\ttotal\\n";\nfor ($i = 0; $i < scalar (@count_motifs); $i++)\n  {\n  my $typ = length ($count_motifs[$i]);\n  print OUT $count_motifs[$i];\n  for ($j = $min_repeats; $j <= $max_repeats; $j++)\n    {\n    if ($j < $typrep{$typ}) {print OUT "\\t-";next};\n    if ($count_motifs[$i]->{$j}) {print OUT "\\t$count_motifs[$i]->{$j}"} else {print OUT "\\t"};\n    };\n  print OUT "\\t$count_motifs{$count_motifs[$i]}\\n";\n  };\nprint OUT "\\n";\n\n#§§§ Frequency of SSRs: summarizing redundant and reverse motifs §§§#\n# Eliminates %count_motifs !\nprint OUT "Frequency of classified repeat types (considering sequence complementary)\\n-------------------------------------------------------------------------\\n\\nRepeats";\nmy (%red_rev,@red_rev); # groups\nfor ($i = 0; $i < scalar (@count_motifs); $i++)\n  {\n  next if ($count_motifs{$count_motifs[$i]} eq \'X\');\n  my (%group,@group,$red_rev); # store redundant/reverse motifs\n  my $reverse_motif = $actual_motif = $actual_motif_a = $count_motifs[$i];\n  $reverse_motif =~ tr/ACGT/TGCA/;\n  $reverse_motif = reverse $reverse_motif;\n  my $reverse_motif_a = $reverse_motif;\n  for ($j = 0; $j < length ($count_motifs[$i]); $j++)\n    {\n    if ($count_motifs{$actual_motif}) {$group{$actual_motif} = "1"; $count_motifs{$actual_motif}=\'X\'};\n    if ($count_motifs{$reverse_motif}) {$group{$reverse_motif} = "1"; $count_motifs{$reverse_motif}=\'X\'};\n    $actual_motif =~ s/(.)(.*)/$2$1/;\n    $reverse_motif =~ s/(.)(.*)/$2$1/;\n    $actual_motif_a = $actual_motif if ($actual_motif lt $actual_motif_a);\n    $reverse_motif_a = $reverse_motif if ($reverse_motif lt $reverse_motif_a)\n    };\n  if ($actual_motif_a lt $reverse_motif_a) {$red_rev = "$actual_motif_a/$reverse_motif_a"}\n  else {$red_rev = "$reverse_motif_a/$actual_motif_a"}; # group name\n  $red_rev{$red_rev}++;\n  @group = keys %group;\n  for ($j = 0; $j < scalar (@group); $j++)\n    {\n    for ($k = $min_repeats; $k <= $max_repeats; $k++)\n      {\n      if ($group[$j]->{$k}) {$red_rev->{"total"} += $group[$j]->{$k};$red_rev->{$k} += $group[$j]->{$k}}\n      }\n    }\n  };\nfor ($i = $min_repeats; $i <= $max_repeats; $i++) {print OUT "\\t$i"};\nprint OUT "\\ttotal\\n";\n@red_rev = sort {length ($a) <=> length ($b) || $a cmp $b } keys %red_rev;\nfor ($i = 0; $i < scalar (@red_rev); $i++)\n  {\n  my $typ = (length ($red_rev[$i])-1)/2;\n  print OUT $red_rev[$i];\n  for ($j = $min_repeats; $j <= $max_repeats; $j++)\n    {\n    if ($j < $typrep{$typ}) {print OUT "\\t-";next};\n    if ($red_rev[$i]->{$j}) {print OUT "\\t",$red_rev[$i]->{$j}}\n    else {print OUT "\\t"}\n    };\n  print OUT "\\t",$red_rev[$i]->{"total"},"\\n";\n  };\n\n''')
    with open('./misa.ini','w') as f:
        f.write('''definition(unit_size,min_repeats):                   1-10 2-6 3-5 4-5 5-5 6-5\ninterruptions(max_difference_between_2_SSRs):        100\n''')



def runCom_Cmd():
    try:
        checkPerl()
        checkCdhit()
        checkMuscle()
        if checkPerl() is True and checkMuscle() is True and checkCdhit() is True:
            showinfo(
                title='info',
                message='Satrt...\nThe run time depends on your file size and computer performance.')

            t1 = time.perf_counter()
            fixGeneId()
            samples = getSamples('.fa')
            creatMisa()
            callMisa(samples)
            getSsrSeq(samples)
            callCdHit()
            reformSSR(samples)
            getReverseSeq()
            #enrichSSR()
            getFinal(enrichSSR())
            t2 = time.perf_counter()

            showinfo(title='info', message='Enrich SSR Done! \n\nTotal time: ' + str(t2 - t1)[:-7] +
                    ' S' '\n\nthe result saved in  Result drictory\n\n')
            infoWindow.delete(0.0, END)
            with open("SSR_Result.txt", "r") as result:
                for line in result:
                    infoWindow.insert('insert', line)
        else:
            raise


    except:
        showinfo(title='Error', message='please check the running environment and Fasta dirctory!')


top = Tk()
top.title('SSR Enricher')
top.geometry('883x565')
style = Style()




VScroll1 = Scrollbar(top, orient='vertical')
VScroll1.place(relx=0.96, rely=0.156, relwidth=0.019, relheight=0.823)

infoWindow = Text(top, )
infoWindow.place(relx=0.263, rely=0.156, relwidth=0.69, relheight=0.823)
infoWindow.insert(
    0.0, "\nSSRENRICHER V1.1\n\nThis software can automatically enrich  polymorphic SSRs.\n\n"
    "You just need to select all of your fasta files from Trinity ( make sure the file suffix is .fasta ),\n\n"
    "And then click the RUN button.\n\n"
    "Meanwhile,  you can make a cup of coffee and wait for the program to finish."
    "\n\n\n\nWarning:You need to install the following software and add to environment variables first："
    "\nBiopython Perl  CD-HIT  muscle\n")
VScroll1.config(command=infoWindow.yview)
infoWindow.config(yscrollcommand=VScroll1.set)


style.configure('TSelectLabel.TLabel', anchor='w')
SelectLabel = Label(top, text='Please slecet', style='TSelectLabel.TLabel')
SelectLabel.place(relx=0.018, rely=0.156, relwidth=0.21, relheight=0.823)


selectComVar = StringVar(value='Select your files')
style.configure('TselectCom.TButton')
selectCom = Button(
    top,
    text='Select',
    textvariable=selectComVar,
    command=selectCom_Cmd,
    style='TselectCom.TButton')
selectCom.place(relx=0.027, rely=0.028, relwidth=0.164, relheight=0.101)


runComVar = StringVar(value='Run')
style.configure('TrunCom.TButton')
runCom = Button(
    top,
    text='Run',
    textvariable=runComVar,
    command=runCom_Cmd,
    style='TrunCom.TButton')
runCom.place(relx=0.381, rely=0.028, relwidth=0.445, relheight=0.101)


top.mainloop()
