import sys
import xml.etree.ElementTree as ET

#NOTE: do we need function to parameterize curve, areas ? -> yes: to implement


def testIfPoint(w):
    w[0]=w[0][0:] #remove first M character
    current=(w[0],w[1])
    print current[0], current[1]   
    if float(w[3])-float(w[0])<=1 or float(w[4])-float(w[2])<=1: #note: si point puis curve n'est pas pris en compte
      return True;
    else:
      return False;

def singularize():
  treeIn = ET.parse('test_medium.svg')
  rootIn=treeIn.getroot()
  print rootIn.attrib
  
  treePoints = ET.parse('points.svg')
  rootPoints=treePoints.getroot()
  
  treeCurves = ET.parse('curves.svg')
  rootCurves=treeCurves.getroot()
  
  treeAreas = ET.parse('area.svg')
  rootAreas =treeAreas.getroot()

  for path in rootIn.iter("path"): #parse all of the path element
    print "parsing the input Image file \n"
    fill = path.attrib['fill'] #irradiance
    d = path.attrib['d'] #distribution
    words=d.split() #list de tous les mots dans d
    
    if words[-1] == "Z": #a mettre dans le fichier area
      print "mot a boucle \n"
      ET.SubElement(rootAreas, 'path', path.attrib)
          
    elif testIfPoint(words)==true: #Point sigularities 
      pointAttribute={"cx":'"'+float(w[0][0:])+'"', "cy":'"'+float(w[1])+'"', "r":"1", "fill":'"'+fill+'"'}   
      ET.SubElement(rootAreas, 'circle', attrib=pointAttribute)
   
      
    else: #curve singularities TO DO
      pass
    
    
    #test if we are in a curve or in an area
    
  treeAreas.write('area.svg')
  treePoints.write('point.svg')
  return;

singularize();
 
 