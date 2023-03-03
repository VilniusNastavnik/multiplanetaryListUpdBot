import re
import sqlite3 as sl
import codecs
from astroquery.simbad import Simbad
import warnings
from datetime import datetime
import time
import logging
import requests

conn = sl.connect('stars.db')
cursor = conn.cursor()

def deg_to_hms(grad,type):

   if(type == "RA"):
     h=float(grad)/15
   else:
     h=float(grad)
	
   m=abs((h-int(h))*60)	
   s=round((m-int(m))*60)
   if(s == 60):
      s = 0
      m +=1
   
   return str(int(h))+"|"+str(int(m))+"|"+str(s)

def hms_to_numb(hms):

   hmsFields = hms.split("|")
   return int(hmsFields[0])*3600 + int(hmsFields[1])*60 + int(hmsFields[2])

def zeroIfEmpty(val):
   """Return 0 if val ='' """
   if(val == ''):
      return 0
   return val

def wikiRightFormat(val,precision):
   """Return '' if val = 0. if precision is specified (>= 0), apply it"""

   if(val == 0 or val == None):
      return ''

   try:
      float(val)
   except: # not convertible to value. Leave it as a string
      return val

   if(precision > 0):
      return str(round(val,precision))
   elif(precision == 0):
      return str(int(round(val)))
   else:
      return str(val)

def getCoordFromSimbad(star):

   if("'s" in star):   #Like Teegarden's
      star = star.replace("'s", "")
   star = re.sub("\s[A-D]$", "", star) #Some binaries (HD 116029 A, HD 177830 A,...) are not found without eliminating the last letter

   cursor.execute("SELECT ra,dec,dist,mag,ids FROM simbad WHERE ids LIKE '%"+star+"%'") 
   rows = cursor.fetchall()
   
   if(len(rows) == 0):
      return None, None, None

   for row in rows:
      if(re.search(star+"[^0-9]",row[4])): #Exclude extra number to avoid confusion between A-123 and A-12
         break

   ra = row[0].split(" ")
   ra[2] = round(float(ra[2]))
   if(ra[2] == 60):
      ra[1] = int(ra[1]) + 1
      ra[2] = 0
   if(ra[1] == 60):
      ra[0] = int(ra[0]) + 1
      ra[1] = 0
   dec = row[1].replace('+', '').split(" ")
   dec[2] = round(float(dec[2]))
   if(dec[2] == 60):
      dec[1] = int(dec[1]) + 1
      dec[2] = 0
   if(dec[1] == 60):
      dec[0] = int(dec[0]) + 1
      dec[1] = 0

   try:
      distance = round(row[2]*3.261563777,1)
   except:
      distance = 0

   return distance, str(int(ra[0]))+"|"+str(int(ra[1]))+"|"+str(ra[2]), str(int(dec[0]))+"|"+str(int(dec[1]))+"|"+str(str(dec[2]))

def getCoordFromSimbadOnline(star):

   name = star
   #stars with peculiar names, not easy to find
   if(star == "1RXS 1609"):
      name = "1RXS J160929.1-210524"
   elif(star == "1SWASP J1407"):
      name = "1SWASP J140747.93-394542.6"
   elif(star == "2M 0103-55 (AB)"):
      name = "SCR J0103-5515"
   elif(star == "Teegarden's"):
      name = "gat 1370"
   elif(star == "Mu Arae"):
      name = "mu Ara" # to avoid confusion with  V* MU Ara
   elif(star == "TOI-4481"):
      name = "GJ 806"
   elif(star == "Nu Ophiuchi"): # to avoid confusion with  V* MU Ara
      name = "HD 163917" 
   elif(star == "PSR 1257 12"):
      name = "PSR B1257+12" 

   Simbad.add_votable_fields('distance')
   result = Simbad.query_object(name)
   distance = 0
   if(result):
      distance = round(float(result['Distance_distance'][0])*3.261563777,1)  # Convert parsec to light years
      raRaw = result['RA'][0].split()
      decRaw = result['DEC'][0].split()

      rh = int(raRaw[0])
      rm = int(raRaw[1])
      rs = round(float(raRaw[2]))
      if(rs == 60):
         rs = 0
         rm += 1 
      if(rm == 60):
         rm = 0
         rh += 1 
    
      dh = int(decRaw[0])
      dm = int(decRaw[1])
      ds = round(float(decRaw[2]))
      if(ds == 60):
         ds = 0
         dm += 1 
      if(dm == 60):
         dm = 0
         dh += 1 

      ra = str(rh)+"|"+str(rm)+"|"+str(rs)
      dec = str(dh)+"|"+str(dm)+"|"+str(ds)
   else:
      #print (star,"not found in Simbad")
      ra = None
      dec = None

   return distance, ra, dec

def getDBRow(name):
  cursor.execute("SELECT name,ra,dec,mag,dist,type,mass,radius,temp,age,metall,planets FROM stars WHERE name = ?",(name, )   )
  row = cursor.fetchone()
  print(row)

def getDataFromExoplanet(exoplanetLocalFile):
   """ Get data from Exoplanet.eu catalog and load them on stars.db """

   if(exoplanetLocalFile == None):
      response = requests.request("GET", "http://exoplanet.eu/catalog/csv/")
      if not response.ok:
         print("Download from http://exoplanet.eu failed")
         exit()
      LinesExo = response.text.splitlines(True)
   else:
      try:
         fileExo = open(exoplanetLocalFile, "r")
      except:
         print("File with data from exoplanet.eu not found")
         exit(0)
   
      next(fileExo) # Skip first line with headers
      LinesExo = fileExo.readlines()
   
   planets = 0  # Variable incremented on every occurence of the same star
   star = ''
   dist=0
   for lineExo in LinesExo:
      if "\"" in lineExo: # There are fields of type "a, b" that are confused with two different fields. This fix it.
         tmpExo = lineExo.replace(', ',';')
      else:
         tmpExo = lineExo
      fieldEx = tmpExo.split(",")

      #read check
      if(len(fieldEx) < 93):
         print("Error in data from exoplanets: expected at least 92 fields, found",len(fieldEx))
         print(tmpExo)
         exit()

      fieldEx[68] = fieldEx[68].rstrip()
      if(fieldEx[68] and fieldEx[68] != "Sun" and re.search(' [a-z]$',fieldEx[0]) is not None): #only valid stars: no the Sun, no brown dwarfs

         if(star and star != fieldEx[68]) : # if it's a new star, save data of the previous

            if(star == "2M0838+15"):
               planets = 0 # Triple brown dwarfs system 
            elif(star == "V470 Cam (AB)"):
               planets = 0 # Orbiting objects are brown dwarfs  

            if(planets > 1):
               #if(star == "55 Cnc"): 
               #   print() 
               sql = '''INSERT INTO stars (name,ra,dec,mag,dist,type,mass,radius,temp,age,metall,planets) VALUES(?,?,?,?,?,?,?,?,?,?,?,?);'''
               try:
                  cursor.execute(sql,(star,ra,dec,mag,dist,type,mass,radius,temp,age,met,planets))
               except sl.Error as err:
                  print("Insert star "+star+"("+ra,dec+") failed:",err)
            planets = 0
            dist = 0

         #save data for next loop           
         ra = deg_to_hms(fieldEx[69],"RA")
         dec = deg_to_hms(fieldEx[70],"DEC")

         if(fieldEx[76]) != '':
            dist = round(float(fieldEx[76])*3.261563777,1)  # Convert parsec to light years
         star=fieldEx[68]
         type=fieldEx[88]
         mag = zeroIfEmpty(fieldEx[71])
         mass = zeroIfEmpty(fieldEx[82])
         radius = zeroIfEmpty(fieldEx[85])
         temp = zeroIfEmpty(fieldEx[92])
         age = zeroIfEmpty(fieldEx[89])
         met = zeroIfEmpty(fieldEx[79])
         planets += 1

   # Last line
   if( planets > 1):
      try:
         cursor.execute(sql,(star,ra,dec,mag,dist,type,mass,radius,temp,age,met,planets))
      except sl.Error as err:
         print(sql,err)

   # Systems listed by NASA and not by Exoplanet.eu
   #try:
   #  cursor.execute(sql,("K2-352","9|21|47","18|28|10",0,0,"",0,0,0,0,0,3))
   #except sl.Error as err:
   #  print(sql,err)

   conn.commit()

def getDataFromSimbad():

   #Simbad.add_votable_fields('distance','flux(V)','measurements','ids')
   Simbad.add_votable_fields('distance','flux(V)','ids')
   
   cursor.execute("SELECT name FROM stars")
   starsRows = cursor.fetchall()
   query = []
   for row in starsRows:
      name = row[0]
      name = re.sub("\s[A-D]$", "", row[0]) #Some binaries (HD 116029 A, HD 177830 A,...) are not found without eliminating the last letter
      query.append(name)

   #stars with peculiar names, not easy to find
   query.append("1RXS J160929.1-210524")
   query.append("1SWASP J140747.93-394542.6")
   query.append("SCR J0103-5515")
   query.append("gat 1370")
   query.append("GJ 806")
   query.append("mu Ara")
   query.append("HD 20794")

   result = Simbad.query_objects(query)
   for row in result:
      sql = '''INSERT INTO simbad (name,ra,dec,mag,dist,ids) VALUES(?,?,?,?,?,?);'''

      distance = row['Distance_distance']
      if(row['Distance_unit'] == 'kpc'):
         distance = distance*1000
      elif(row['Distance_unit'] == 'Mpc'):
         distance = distance*1000000
      try:
         cursor.execute(sql,(row['MAIN_ID'],row['RA'],row['DEC'],row['FLUX_V'],distance,row['IDS']))
      except sl.Error as err:
         if("UNIQUE constraint failed" not in err.args[0]): #Ignore duplicated
            print("Insert star "+row['MAIN_ID']+"("+row['RA'],row['DEC']+") failed:",err)
   
   conn.commit()

   #update stars
   sql = "UPDATE stars SET dist=COALESCE(NULLIF(dist,''),?),ra=?,dec=? WHERE name=?;"
   for row in starsRows:
      dist, ra, dec = getCoordFromSimbad(row[0])
      if(ra == None): # if not found in Simbad local try again online
         dist, ra, dec = getCoordFromSimbadOnline(row[0])

      if(ra != None):
         try:
            cursor.execute(sql,(dist,ra,dec,row[0]))
         except sl.Error as err:
            print("Update 'stars' table with Simbad data for",row[0],"failed:",err)

   conn.commit()

def getDataFromNASA(nasaLocalFile):
   """ Get data from NASA's catalog and add them on stars.db """

   if(nasaLocalFile == None):
      response = requests.request("GET", "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+hostname,ra,dec,sy_vmag,sy_dist,st_spectype,st_mass,st_rad,st_teff,st_age,st_met+from+pscomppars&format=csv")
      if not response.ok:
         print("Download from https://exoplanetarchive.ipac.caltech.edu failed")
         exit()
      #linesNASAtmp = response.text.splitlines(True)
      #linesNASA = (list(set(linesNASAtmp)))
      linesNASA = (list(set(response.text.splitlines(True)))) #split resulting text in lines and eliminate duplicates
   else:
      try:
         fileNASA = open(nasaLocalFile, "r")
      except:
         print("File with data from exoplanetarchive.ipac.caltech.edu not found")
         exit(0)
   
      next(fileNASA) # Skip first line with headers
      linesNASA = fileNASA.readlines()

   for lineNASA in linesNASA:   
      ts = datetime.timestamp(datetime.now())
      fieldNASA = lineNASA.strip("'").rstrip('\n').split(",")
    
      name = fieldNASA[0].strip('"')
      if(name == "hostname"): #Skip line with headers
         continue

      dist = 0
      if(fieldNASA[4] != ''):
         try:
            dist = round(float(fieldNASA[4])*3.261563777,1)  # Convert parsec to light years
         except:
            print(fieldNASA[4])
            print()

      cursor.execute("SELECT * FROM stars WHERE name = ?",(name, ) )
      rowForName = cursor.fetchone()
      if(rowForName is None):# name not found. Try to update by coordinates
         distSimbad, ra, dec = getCoordFromSimbad(name)
         sql = "UPDATE stars SET mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,''),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,''),?),radius=COALESCE(NULLIF(radius,''),?),temp=COALESCE(NULLIF(temp,''),?),age=COALESCE(NULLIF(age,''),?),metall=COALESCE(NULLIF(metall,''),?)  WHERE ra=? AND dec=?;"
         try:
            cursor.execute(sql,(fieldNASA[3],dist,fieldNASA[5],fieldNASA[6],fieldNASA[7],fieldNASA[8],fieldNASA[9],fieldNASA[10],ra, dec))
         except sl.Error as err:
            print("Update NASA data to",name,"failed:",err)
      else: #udate by name
         sql = "UPDATE stars SET mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,''),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,''),?),radius=COALESCE(NULLIF(radius,''),?),temp=COALESCE(NULLIF(temp,''),?),age=COALESCE(NULLIF(age,''),?),metall=COALESCE(NULLIF(metall,''),?)  WHERE name=?;"
         try:
            cursor.execute(sql,(fieldNASA[3],dist,fieldNASA[5],fieldNASA[6],fieldNASA[7],fieldNASA[8],fieldNASA[9],fieldNASA[10],name))
         except sl.Error as err:
            print("Update 'stars' tables with NASA data for ",name,"failed:",err)

   conn.commit()

def getDataFromWikipedia(wikiLocalFile):
   """ Get data from Wikipedia's page and add them on stars.db """

   if(wikiLocalFile == None):
      response = requests.request("GET", "https://it.wikipedia.org/w/index.php?action=raw&title=Sistemi_multiplanetari")
      if not response.ok:
         print("Download from https://it.wikipedia.org/w/index.php?action=raw&title=Sistemi_multiplanetari failed")
         exit()
      fileWikiRaw = response.text.splitlines(True)
   #else:
   #   try:
   #      fileNASA = open(wikiLocalFile, "r")
   #   except:
   #      print("File with data from exoplanetarchive.ipac.caltech.edu not found")
   #      exit(0)

   #with codecs.open(wikiLocalFile, encoding='utf-8') as fileWikiRaw:
      for lineWikiRaw in fileWikiRaw:
         if "Stella" in lineWikiRaw:  # put wiki data in fieldWR 
            tmpWiki = lineWikiRaw.replace('||',';').replace(',','.').replace('−','-')[1:].rstrip('}}\n')  #Template sep || -> ; and dec sep , -> .
            fieldWR = tmpWiki.split(";")

            #Split values from names
            valsWiki = []
            for i in range(11):
               tmp = fieldWR[i].split("=")[1]
               if(tmp != '-' and tmp != ''):
                  valsWiki.append(tmp)
               else:
                  valsWiki.append(None)

            wikiName = valsWiki[0].replace('[[','').replace(']]','').encode("utf-8").decode()  #get star name from wiki page
            tmpName = wikiName.split("|")
            name = tmpName[0] # Appearing name, not internal link

            raWf = valsWiki[1].split("|")
            ma = int(raWf[2])
            sa = round(float(raWf[3][:-2]))
            if(sa == 60):
               sa = 0
               ma +=1
            raW = str(int(raWf[1])) +"|"+ str(ma) +"|"+ str(sa)      #get RA from wiki page
            
            decWf = valsWiki[2].split("|")
            md = int(decWf[2])
            sd = round(float(decWf[3][:-2]))
            if(sd == 60):
               sd = 0
               md +=1
            decW = str(int(decWf[1])) +"|"+ str(md) +"|"+ str(sd)    #get DEC from wiki page
            
            # Search coordinates in 'stars' by Simbad(name) (distance ignored ' cause alredy in 'stars')
            distSimbad, ra, dec = getCoordFromSimbad(name)
            if(ra == None): # if not found in Simbad local try again online
               distSimbad, ra, dec = getCoordFromSimbadOnline(name)

            if(ra == None): # if not found anyway, use wiki data
               ra = raW
               dec = decW

            #search by coordinates in 'stars' db
            cursor.execute("SELECT name,mag,dist,type,mass,radius,temp,age,metall FROM stars WHERE ra=? AND dec=?",(ra,dec))
            row = cursor.fetchone()
            if(row is None): # or search by name
               try:
                  cursor.execute("SELECT name,mag,dist,type,mass,radius,temp,age,metall FROM stars WHERE name='"+name+"'")
                  row = cursor.fetchone()
               except sl.Error as err:
                  print(name,err)

            # UPDATE WHERE not empty or not too different
            if(row is None):  #It's in wiki, not in 'stars' db. Maybe ther is a problem
               print("Wiki star",name,"with coordinates RA:",raW,"DEC:",decW,"not in catalog exoplanet/Simbad") 
            else:
               # Use wiki name and data if not null
               rowList = list(row)

               distStars = rowList[2]  #dist in 'stars'
               dist = valsWiki[4]      #dist in Wiki
               if(distStars != "" and distStars != None and distStars > 0):
                  if(abs ((float(valsWiki[4]) - float(distStars))) > 1000):  # Too different. Maybe there is a problem
                     print(name,"diff distance exo-wiki ("+str(distStars)+"-"+str(valsWiki[4])+") > 100. Keep Wiki")
                  else:
                     dist = distStars
               
               sqlUpd = "UPDATE stars SET name=?,mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,0),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,0),?),radius=COALESCE(NULLIF(radius,0),?),temp=COALESCE(NULLIF(temp,0),?),age=COALESCE(NULLIF(age,0),?),metall=COALESCE(NULLIF(metall,0),?) WHERE ra=? AND dec=?"   
               cursor.execute(sqlUpd, (wikiName,valsWiki[3],dist,valsWiki[5],valsWiki[6],valsWiki[7],valsWiki[8],valsWiki[9],valsWiki[10],ra,dec) )

   conn.commit()

def generateWikitable(tableOutFile):
   """ Generate wikitables from on stars.db """
   tableWiki = open(tableOutFile, 'w', encoding="utf-8")

   cursor.execute("SELECT COUNT(*), planets FROM stars GROUP BY planets ORDER BY planets ASC")
   rows = cursor.fetchall()

   p = [0,0,0,0,0,0,0,0,0]
   for i in range(len(rows)):
      p[i] = rows[i][0]

   line = "{{Progetto sistemi multiplanetari|con2pianeti="+str(p[0])+"|con3pianeti="+str(p[1])+"|con4pianeti="+str(p[2])+"|con5pianeti="+str(p[3])+"|con6pianeti="+str(p[4])+"|con7pianeti="+str(p[5])+"|con8pianeti="+str(p[6])+"}}"
   tableWiki.write("%s\n" % line)

   tableWiki.write("%s\n" % "")

   line = "<noinclude>{{Stelle con pianeti extrasolari confermati/Top}}</noinclude>"
   tableWiki.write("%s\n" % line)

   cursor.execute("SELECT * FROM stars ORDER BY name ASC")
   rows = cursor.fetchall()
   for row in rows:
      line = "{{Stelle con pianeti extrasolari confermati"
      tableWiki.write("%s\n" % line)
      raIn = row[1].split("|")
      decIn = row[2].split("|")
      raOut = raIn[0].zfill(2)+"|"+raIn[1].zfill(2)+"|"+raIn[2].zfill(2)
      if(int(decIn[0]) < 0):
         tmp = "-"+str(-int(decIn[0])).zfill(2)
      else:
         tmp = decIn[0].zfill(2)
      decOut = tmp+"|"+decIn[1].zfill(2)+"|"+decIn[2].zfill(2)

      mag = wikiRightFormat(row[3],2)
      dist = wikiRightFormat(row[4],-1)
      tipo = wikiRightFormat(row[5],-1)
      massa = wikiRightFormat(row[6],2)
      try:
         raggio = wikiRightFormat(row[7],2)
      except:
         print(row[0])
      temper = wikiRightFormat(row[8],0)
      eta = wikiRightFormat(row[9],-1)
      metall = wikiRightFormat(row[10],2)

      tmp = "|Stella=[["+row[0]+"]]||Ascensione retta={{RA|"+raOut+"}}||Declinazione={{DEC|"+decOut+"}}||Magnitudine apparente="+mag+"||Distanza="+dist+"||Tipo spettrale="+tipo+"||Massa="+massa+"||Raggio="+raggio+"||Temperatura="+temper+"||Età="+eta+"||Metallicità="+metall+"||Pianeti="+str(row[11])+"}}"
      
      tableWiki.write("%s\n" % tmp.replace('.',','))

   line = "<noinclude>{{Stelle con pianeti extrasolari confermati/Bottom}}</noinclude>"
   tableWiki.write("%s\n" % line)

   conn.close()

def main():

   warnings.simplefilter('ignore', UserWarning)

   cursor.execute("DROP TABLE IF EXISTS stars")
   with conn:
      conn.execute("""
         CREATE TABLE stars (
         name TEXT,
         ra TEXT,
         dec TEXT,
         mag REAL DEFAULT 0.0,
         dist REAL DEFAULT 0.0,
         type TEXT,
         mass REAL DEFAULT 0.0,
         radius REAL DEFAULT 0.0,
         temp REAL DEFAULT 0.0,
         age REAL DEFAULT 0.0,
         metall REAL DEFAULT 0.0,
         planets INTEGER
      );
      """)

   cursor.execute("DROP TABLE IF EXISTS simbad")
   with conn:
      conn.execute("""
         CREATE TABLE simbad (
         name TEXT,
         ra TEXT,
         dec TEXT,
         mag REAL,
         dist REAL,
         type TEXT,
         mass REAL,
         radius REAL,
         temp REAL,
         age REAL,
         metall REAL,
         ids TEXT,
         PRIMARY KEY (ra, dec)
      );
      """)

   print("Retrieving data from Exoplanet ...")
   getDataFromExoplanet(None)    
   #getDataFromExoplanet("exoplanet.out")
   print("Retrieving data from Simbad ...")
   getDataFromSimbad()
   print("Retrieving data from NASA ...")
   getDataFromNASA(None)                     
   #getDataFromNASA("NASA.out")
   print("Retrieving data from Wikipedia ...")
   getDataFromWikipedia(None)
   #getDataFromWikipedia("wiki.out")
   print("Generating 'tabella.wiki' ...")
   generateWikitable("tabella.wiki")
   print("Done. Copy the content of the file 'tabella.wiki' in the correct place inside https://it.wikipedia.org/wiki/Sistemi_multiplanetari. Check the result before publishing!")

if __name__ == "__main__":
   main()