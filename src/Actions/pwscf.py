import sys
import os
import commands

def create(infile_pre, infile_post, coordfile):
  filename = 'pwscf_pimcpp.in'
  runfile = open(filename,'w')
  i = open(infile_pre,'r')
  line = i.readline()
  while(line != ''):
    runfile.write(line)
    line = i.readline()
  i.close()
  X = open(coordfile,'r')
  line = X.readline()
  while(line != ''):
    runfile.write(line)
    line = X.readline()
  X.close()
  #for x in range(0,len(coords)):
  #    runfile.write(coords[x])
  #    #runfile.write(coords[x] +'\n')
  i = open(infile_post,'r')
  line = i.readline()
  while(line != ''):
    runfile.write(line)
    line = i.readline()
  i.close()
  print "Created pwscf input",filename
  return filename

def run(filename):
  print "Running pwscf..."
  os.system('mpirun -np 4 pw.x < ' + filename + ' > pw.out')
  x = commands.getoutput('grep energy pw.out')
  print "pwscf finished."
  return x

def processString(line):
  entry = []
  string = ''
  i = 0
  while i < len(line):
    if (line[i] == ' ' or line[i] == ',' or line[i] == '\t' or line[i] == '\n'):  #Edit for the delimiters you want to recognize
      if (line[i-1] != ' ' and i>0):
        entry.append(string)
        string = ''
    else:
      string += line[i]
    i += 1
  return entry

def parse(output):
  p = processString(output)
  e = []
  index = 0
  while(index<len(p)):
    if(p[index] == 'total'):     
      index += 3
      if(p[index] != 'the'):
        e.append(float(p[index]))
        #print p[index]
    index+=1
  #print "collected energies",e
  return(e[len(e)-1])

def execute(infile_pre, infile_post, coords, energy_out):
  f = create(infile_pre, infile_post, coords)
  output = run(f)
  converged_E = parse(output)
  out=open(energy_out,'w')
  out.write(str(converged_E))
  out.close()
