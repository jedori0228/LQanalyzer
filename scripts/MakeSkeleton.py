import os,sys

NAME = sys.argv[1]

LQDIR = os.environ['LQANALYZER_DIR']

SkeletonDir = LQDIR+'/JskimData/AnalyzerSkeleton/'

lines_h = open(SkeletonDir+'AnalyzerSkeleton.h').readlines()
lines_cc = open(SkeletonDir+'AnalyzerSkeleton.cc').readlines()

AnalyzerDir = LQDIR+'/LQAnalysis/Analyzers/'

new_h = open(AnalyzerDir+'include/'+NAME+'.h','w')
for line in lines_h:
  if 'AnalyzerSkeleton' in line:
    line = line.replace('AnalyzerSkeleton', NAME)
  new_h.write(line)
new_h.close()

new_cc = open(AnalyzerDir+'src/'+NAME+'.cc','w')
for line in lines_cc:
  if 'AnalyzerSkeleton' in line:
    line = line.replace('AnalyzerSkeleton', NAME)
  new_cc.write(line)
new_cc.close()







