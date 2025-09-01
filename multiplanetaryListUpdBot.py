import re
import sqlite3 as sl
from astroquery.simbad import Simbad
import warnings
from astroquery.exceptions import AstropyWarning
import logging
import requests
#import codecs
#import time
#from datetime import datetime

sqliteConn = sl.connect('stars.db')
sqliteCursor = sqliteConn.cursor()
brownDwarfMassLimit = 13.0

peculiar_star_names = {  #stars with peculiar names for Simbad, not easy to find
        "1RXS 1609": "1RXS J160929.1-210524",
        "1SWASP J1407": "1SWASP J140747.93-394542.6",
        "2M 0103-55 (AB)": "SCR J0103-5515",
        "Teegarden's": "gat 1370",
        "Mu Arae": "mu Ara",  # to avoid confusion with V* MU Ara
        "Nu Ophiuchi": "HD 163917",  # to avoid confusion with V* NU Oph
        "82 Eri": "HD 20794"
}

def deg_to_hms(grad,cooType):

    logging.debug(f"deg_to_hms: {grad}, {cooType}")
    
    if(cooType == "RA"):
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

def getCoordFromSimbadLocalTable(star):

    logging.debug(f"getCoordFromSimbad: {star}")

    sqliteCursor.execute("SELECT ra, dec, dist, mag, ids FROM simbad WHERE REPLACE(ids, '  ', ' ') LIKE '%' || REPLACE(?, '  ', ' ') || '%'", (star,))
    rows = sqliteCursor.fetchall()
   
    if(len(rows) == 0):
        return None, None, 0, 0, None
    if(len(rows) == 1):
        return rows[0][0], rows[0][1], rows[0][2], rows[0][3], rows[0][4]

    # In case of multiple records
    count = 0
    pattern = re.escape(star) + r'\||' + re.escape(star) + r'$'
    for row in rows:
        ids = re.sub(r'\s+', ' ', row[4])  # Remove multiple spaces
        if re.search(pattern, ids, re.IGNORECASE):
            rec = row
            count += 1

    if count != 1:
        logging.warning(f"getCoordFromSimbadLocalTable: found {count} data record for {star}.")
        return None, None, 0, 0, None
    else:
        return rec[0], rec[1], rec[2], rec[3], rec[4]

def query_simbad(name):
    result = Simbad.query_object(name)
    if result: # if found use only the first result [0]
        ra = deg_to_hms(result['ra'][0], "RA")
        dec = deg_to_hms(result['dec'][0], "DEC")
        distance = pc2LigthYear(result['mesdistance.dist'][0], result['mesdistance.unit'][0])
        mag = result['V'][0]
        ids = re.sub(r'\|\*\s{1,2}', '|', result['ids'][0])  #drop leading asterisks and spaces
        ids = re.sub(r'NAME\s+', '', ids)  # Remove "NAME " from the name
        return ra, dec, distance, mag, ids
    return None, None, 0, 0, None

def pc2LigthYear(distance, unit):
    distance = round(float(distance) * 3.261563777, 1)  # Convert parsecs to light-years
    if 'kpc' in unit:
        distance = distance*1000
    elif 'Mpc' in unit:
        distance = distance*1000000

    return distance
    
def getCoordFromSimbadOnline(star):

    logging.debug(f"getCoordFromSimbadOnline: {star}")
    
    name = star
    if star in peculiar_star_names:
        name = peculiar_star_names.get(star)

    Simbad.add_votable_fields('mesdistance','V','ids')    

    if re.search(r"\s[A-D]$", star):
        name = re.sub(r"\s[A-D]$", "", star)  # Remove last letter for multiple systems, sometimes cataloged as single

    if re.search(r"\s\(AB\)$", star):
        name = re.sub(r"\s\(AB\)$", "", star)  # Remove characters indicating a double star, sometimes cataloged as single

    response = query_simbad(star)
    if response[0]:
        return response
    else:
        response = query_simbad(name)
        if response[0]:
            return response
        else:
            # If no results found
            logging.debug(f"getCoordFromSimbadOnline: no results found for {star}")
            return None, None, 0, 0, None

def getDBRow(name):
    sqliteCursor.execute("SELECT name,ra,dec,mag,dist,type,mass,radius,temp,age,metall,planets FROM stars WHERE name = ?",(name, )   )
    row = sqliteCursor.fetchone()
    print(row)

def getDataFromExoplanet(exoplanetLocalFile):
    """ Get data from Exoplanet.eu catalog and load them on stars.db """
    
    logging.info(f"getDataFromExoplanet")

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
   
    planets = 0  # Variable incremented on every occurrence of the same star
    star = ''
    ra=0;dec=0;mag=0;dist=0;spec_type="";mass=0;radius=0;temp=0;age=0;met=0;altNames=""
    for lineExo in LinesExo:
        
        #if "Kepler-451" in lineExo:
        #    print() 
                
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

        fieldEx[68] = fieldEx[68].rstrip()  # Mother star        
        if(fieldEx[68] and fieldEx[68] != "Sun" and re.search(r'[A-Z\)\s][a-z]$',fieldEx[0]) is not None): #only valid stars: no the Sun. It means 'Name b' or 'Name AB)b' (at least one planet)
               
            planet_mass = float(zeroIfEmpty(fieldEx[2]))
            if planet_mass > brownDwarfMassLimit:  # no brown dwarfs
                logging.info(f"{fieldEx[0]} mass is {planet_mass}: excluded as probable brown dwarf")
            else:
                if(star and star != fieldEx[68]) : # if it's a new star, save data of the previous one
                        
                    if(star == "HS 0705+6700"):
                        planets = 0 # Possible brown dwarfs system 
    
                    if(planets > 1):
                        sql = '''INSERT INTO stars (name,ra,dec,mag,dist,type,mass,radius,temp,age,metall,planets,altNames) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?);'''
                        try:
                            sqliteCursor.execute(sql,(star,ra,dec,mag,dist,spec_type,mass,radius,temp,age,met,planets,altNames))
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
                spec_type=fieldEx[88]
                mag = zeroIfEmpty(fieldEx[71])
                mass = zeroIfEmpty(fieldEx[82])
                radius = zeroIfEmpty(fieldEx[85])
                temp = zeroIfEmpty(fieldEx[92])
                age = zeroIfEmpty(fieldEx[89])
                met = zeroIfEmpty(fieldEx[79])                
                altNames = fieldEx[97]
                planets += 1

    # Last line
    if( planets > 1):
        try:
            sqliteCursor.execute(sql,(star,ra,dec,mag,dist,spec_type,mass,radius,temp,age,met,planets,altNames))
        except sl.Error as err:
            print(sql,err)

    # Systems listed by NASA and not by Exoplanet.eu
    #try:
    #  sqliteCursor.execute(sql,("K2-352","9|21|47","18|28|10",0,0,"",0,0,0,0,0,3))
    #except sl.Error as err:
    #  print(sql,err)

    sqliteConn.commit()

def getDataFromSimbadSite():
    """
    Fetches star data from the Simbad astronomical database and updates the local database for coordinates distance (in pc) and magnitude.
    If a star is not found with the first general query, it attempts to retrieve more precise data online.

    """

    logging.info(f"getDataFromSimbadSite") 
   
    #Get data already downloaded from Exoplanets.eu
    sqliteCursor.execute("SELECT name FROM stars")
    starsRows = sqliteCursor.fetchall()
    query = []
    for row in starsRows:
        name = row[0]
        query.append(name)

    query.extend(peculiar_star_names.values())

    Simbad.add_votable_fields('mesdistance','V','ids')  
    result = Simbad.query_objects(query)

    #put Simbad results in local db (table 'simbad')
    for row in result:
        sql = '''INSERT INTO simbad (name,ra,dec,mag,dist,ids) VALUES(?,?,?,?,?,?);'''

        if(row['main_id'] != ''): # if name is empty, it's a wrong record
            name = re.sub(r'^\*\s{1,2}', '', row['main_id']) #drop leading asterisks and spaces
            name = re.sub(r'^NAME\s+', '', name)  # Remove "NAME " from the beginning of the name

            ra = deg_to_hms(row['ra'], "RA")
            dec = deg_to_hms(row['dec'], "DEC")
            distance = pc2LigthYear(row['mesdistance.dist'], row['mesdistance.unit'])
            mag = result['V'][0]
            if not mag:
                mag = 0
            ids = re.sub(r'\|\*\s{1,2}', '|', result['ids'][0])  #drop leading asterisks and spaces
            ids = re.sub(r'NAME\s+', '', ids)  # Remove "NAME " from the name
            ids = re.sub(r'\|\*\s{1,2}', '|', row['ids'])  #drop leading asterisks and spaces
            ids = re.sub(r'NAME\s+', '', ids)  # Remove "NAME " from the name
            ids = name + "|" + ids

            try:
                sqliteCursor.execute(sql,(name,ra,dec,mag,distance,ids))
            except sl.Error as err:
                if("UNIQUE constraint failed" not in err.args[0]): #Ignore duplicated
                    print("Insert star "+name+"("+ra,dec+") failed:",err)
    
    sqliteConn.commit()

    # double check 'simbad' table
    for row in starsRows:
        #Search in 'simbad' table if star is correctly present
        ra, dec, dist, mag, ids = getCoordFromSimbadLocalTable(row[0])
        if(ra == None): # if not found in Simbad local DB try again online. It's more efficient in finding names

            ra, dec, dist, mag, ids = getCoordFromSimbadOnline(row[0])

            if(ra != None): #if found, check if already present with different ids and update 'simbad' table ids. Otherwise insert new row
                try:
                    # Check for star with this new found ids (synonyms)
                    sqliteCursor.execute("SELECT name FROM simbad WHERE ids LIKE ?", ('%' + ids + '%',))
                    synonymIds = sqliteCursor.fetchone()
                    ids = row[0] + "|" + ids
                    update = False
                    if synonymIds: #if found, add the name as new id
                        sqliteCursor.execute("UPDATE simbad SET ids=? WHERE name=?;",(ids,synonymIds[0]))
                        logging.debug(f"Star {synonymIds[0]} data updated with data from Simbad")
                        update = True
                    else: # otherwise check for star with this new found coordinats (synonyms)
                        sqliteCursor.execute("SELECT name, ra, dec FROM simbad WHERE ra=? AND dec=?", (ra,dec))
                        synonymCoord = sqliteCursor.fetchone()
                        if synonymCoord:
                            sqliteCursor.execute("UPDATE simbad SET ids=? WHERE ra=? AND dec=?",(ids,ra,dec))
                            logging.debug(f"Star of coordinates {ra} and {dec} updated with data from Simbad")
                            update = True
                    if not update: # otherwise insert new row
                        if not mag:
                            mag = 0
                        sqliteCursor.execute("INSERT INTO simbad (name,ra,dec,mag,dist,ids) VALUES(?,?,?,?,?,?);",
                                (row[0],ra,dec,mag,dist,ids))
                        logging.debug(f"Star {row[0]} inserted in simbad table")
                except sl.Error as err:
                    print("Update/insert in 'simbad' table with online Simbad data for",row[0],"failed:",err)

        sqliteConn.commit()

    # Update 'star' table with data from 'simbad' table
    for row in starsRows:
        ra, dec, dist, mag, ids = getCoordFromSimbadLocalTable(row[0])
        if(ra != None):
            # If found, update 'star' table
            sqlStars = """
                UPDATE stars
                SET dist = CASE WHEN dist = 0 OR dist IS NULL OR dist = '' THEN ? ELSE dist END,
                    mag = CASE WHEN mag = 0 OR mag IS NULL OR mag = '' THEN ? ELSE mag END,
                    ra = ?,
                    dec = ?
                WHERE name = ?;
            """
            try:
                sqliteCursor.execute(sqlStars,(dist,mag,ra,dec,row[0]))
                logging.debug(f"Star {row[0]} data updated with data from Simbad")
            except sl.Error as err:
                print("Update 'stars' table with Simbad data for",row[0],"failed:",err)
        else:
            logging.info(f"Star {row[0]} not present in Simbad")

    sqliteConn.commit()

def getDataFromNASA(nasaLocalFile):
    """ Get data from NASA's catalog and add them on stars.db """

    logging.info(f"getDataFromNASA")
    
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
        #ts = datetime.timestamp(datetime.now())
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

        sqliteCursor.execute("SELECT * FROM stars WHERE name = ?",(name, ) )
        rowForName = sqliteCursor.fetchone()
        if(rowForName is None):# name not found. Try to update by coordinates
            ra, dec, distSimbad, mag, ids = getCoordFromSimbadLocalTable(name)
            sql = "UPDATE stars SET mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,''),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,''),?),radius=COALESCE(NULLIF(radius,''),?),temp=COALESCE(NULLIF(temp,''),?),age=COALESCE(NULLIF(age,''),?),metall=COALESCE(NULLIF(metall,''),?)  WHERE ra=? AND dec=?;"
            try:
                sqliteCursor.execute(sql,(fieldNASA[3],dist,fieldNASA[5],fieldNASA[6],fieldNASA[7],fieldNASA[8],fieldNASA[9],fieldNASA[10],ra, dec))
                logging.debug(f"Star of coordinates RA:{ra},DEC:{dec} updated with data from NASA")
            except sl.Error as err:
                print("Update NASA data to",name,"failed:",err)
                logging.error(f"Update NASA data to {name} failed: {err}")

        else: #udate by name
            sql = "UPDATE stars SET mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,''),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,''),?),radius=COALESCE(NULLIF(radius,''),?),temp=COALESCE(NULLIF(temp,''),?),age=COALESCE(NULLIF(age,''),?),metall=COALESCE(NULLIF(metall,''),?)  WHERE name=?;"
            try:
                sqliteCursor.execute(sql,(fieldNASA[3],dist,fieldNASA[5],fieldNASA[6],fieldNASA[7],fieldNASA[8],fieldNASA[9],fieldNASA[10],name))
                logging.debug(f"Star {name} updated with data from NASA")
            except sl.Error as err:
                print("Update 'stars' tables with NASA data for ",name,"failed:",err)
                logging.error(f"Update 'stars' tables with NASA data for {name} failed: {err}")

    sqliteConn.commit()

def getDataFromWikipedia(wikiLocalFile):
    """ Get data from Wikipedia's page and add them on stars.db """

    logging.debug(f"getDataFromWikipedia")
    
    if(wikiLocalFile == None):
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        response = requests.get("https://it.wikipedia.org/w/index.php?action=raw&title=Sistemi_multiplanetari",headers)
        if not response.ok:
            print("Download from https://it.wikipedia.org/w/index.php?action=raw&title=Sistemi_multiplanetari failed. Program terminated.")
            logging.error(f"Failed to retrieve data from Wikipedia: {response}")
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
                #if "HD 219134" in valsWiki[0]:
                #    print(valsWiki[0])
                wikiName = valsWiki[0].replace('[[','').replace(']]','').encode("utf-8").decode()  #get star name from wiki page
                tmpName = wikiName.split("|")
                name = tmpName[0] # Appearing name, not internal link
                
                ref_index = name.find("<ref>")
                if ref_index != -1: #if the name contains a reference, exclude it
                    name = name[0:ref_index]
                    
                if "<ref>" in lineWikiRaw:
                    logging.warning(f"A reference is present in {name}. Check that it's correctly reported.")
                    print("A reference is present in "+name+". Check that it's correctly reported.")

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
            
                # Search coordinates in 'stars' by Simbad(name) (distance ignored 'cause already in 'stars')
                ra, dec, distSimbad, mag, ids = getCoordFromSimbadLocalTable(name)
                """if(ra == None): # if not found in Simbad local try again online
                    distSimbad, ra, dec = getCoordFromSimbadOnline(name)
                    #print(distSimbad, ra, dec)"""

                if(ra == None): # if not found anyway, use wiki data
                    ra = raW
                    dec = decW

                #search by coordinates in 'stars' db
                sqliteCursor.execute("SELECT name,mag,dist,type,mass,radius,temp,age,metall FROM stars WHERE ra=? AND dec=?",(ra,dec))
                row = sqliteCursor.fetchone()
                if(row is None): # or search by name
                    try:
                        name2Search = name
                        if name.startswith("Gliese"):
                            name2Search = name.replace("Gliese", "GJ", 1)  #Gliese stars appear as GJ in exoplanet.eu
                        #elif name.endswith(" [A-D]"):
                        #    name2Search = name.replace(" [A-D]", "", 1)  #
                        sqliteCursor.execute("SELECT name,mag,dist,type,mass,radius,temp,age,metall FROM stars WHERE name LIKE'"+name2Search+"%'")
                        row = sqliteCursor.fetchone()
                    except sl.Error as err:
                        print(name,err)

                # UPDATE WHERE not empty or not too different
                if(row is None):  #It's in wiki, not in 'stars' db. Maybe ther is a problem
                    logging.info(f"Star {name} with coordinates RA:{raW}, DEC:{decW} is in current wiki, but not valid (not existent or less then 2 valid planets)") 
                    print("Star",name,"with coordinates RA:",raW,"DEC:",decW," is in current wiki, but not valid (not existent or less then 2 valid planets)") 
                else:
                    # Use wiki name and data if not null
                    rowList = list(row)

                    distStars = rowList[2]  #dist in 'stars'
                    dist = valsWiki[4]      #dist in Wiki
                    if(distStars != "" and distStars != None and distStars > 0):
                        if(abs ((float(valsWiki[4]) - float(distStars))) > 1000):  # Too different. Maybe there is a problem
                            print(name,"difference in distance between Exoplanet and Wikipedia is "+str(distStars)+"-"+str(valsWiki[4])+" > 100 ly. Keeping Wikipedia datum")
                            logging.warning(f"difference in distance between Exoplanet and Wikipedia for star {name} is "+str(distStars)+"-"+str(valsWiki[4])+" > 100 ly. Keeping Wikipedia datum")
                        else:
                            dist = distStars
               
                    sqlUpd = "UPDATE stars SET name=?,mag=COALESCE(NULLIF(mag,0),?),dist=COALESCE(NULLIF(dist,0),?),type=COALESCE(NULLIF(type,''),?),mass=COALESCE(NULLIF(mass,0),?),radius=COALESCE(NULLIF(radius,0),?),temp=COALESCE(NULLIF(temp,0),?),age=COALESCE(NULLIF(age,0),?),metall=COALESCE(NULLIF(metall,0),?) WHERE ra=? AND dec=?"   
                    sqliteCursor.execute(sqlUpd, (wikiName,valsWiki[3],dist,valsWiki[5],valsWiki[6],valsWiki[7],valsWiki[8],valsWiki[9],valsWiki[10],ra,dec) )

    sqliteConn.commit()

def generateWikitable(tableOutFile):
    """ Generate wikitables from on stars.db """
    
    logging.debug(f"generateWikitable")
    
    tableWiki = open(tableOutFile, 'w', encoding="utf-8")

    sqliteCursor.execute("SELECT COUNT(*), planets FROM stars GROUP BY planets ORDER BY planets ASC")
    rows = sqliteCursor.fetchall()

    p = [0,0,0,0,0,0,0,0,0]
    for i in range(len(rows)):
        p[i] = rows[i][0]

    line = "{{Progetto sistemi multiplanetari|con2pianeti="+str(p[0])+"|con3pianeti="+str(p[1])+"|con4pianeti="+str(p[2])+"|con5pianeti="+str(p[3])+"|con6pianeti="+str(p[4])+"|con7pianeti="+str(p[5])+"|con8pianeti="+str(p[6])+"}}"
    tableWiki.write("%s\n" % line)

    tableWiki.write("%s\n" % "")

    line = "<noinclude>{{Stelle con pianeti extrasolari confermati/Top}}</noinclude>"
    tableWiki.write("%s\n" % line)

    sqliteCursor.execute("SELECT * FROM stars ORDER BY name ASC")
    rows = sqliteCursor.fetchall()
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

        ref_index = row[0].find("<ref>")
        if ref_index != -1: #if the name contains a reference put it outside [[name]]
            name = "[["+row[0][0:ref_index]+"]]"+row[0][ref_index:]
        else:
            name = "[["+row[0]+"]]"
        tmp = "|Stella="+name+"||Ascensione retta={{RA|"+raOut+"}}||Declinazione={{DEC|"+decOut+"}}||Magnitudine apparente="+mag+"||Distanza="+dist+"||Tipo spettrale="+tipo+"||Massa="+massa+"||Raggio="+raggio+"||Temperatura="+temper+"||Età="+eta+"||Metallicità="+metall+"||Pianeti="+str(row[11])+"}}"
      
        tableWiki.write("%s\n" % tmp.replace('.',','))

    line = "<noinclude>{{Stelle con pianeti extrasolari confermati/Bottom}}</noinclude>"
    tableWiki.write("%s\n" % line)

    sqliteConn.close()

def main():

    warnings.simplefilter('ignore', UserWarning)
    
    # Set up logging configuration
    logging.basicConfig(
        level=logging.INFO,  # Set the minimum log level to DEBUG
        format='%(asctime)s - %(levelname)s - %(message)s',  # Log format
        handlers=[
            logging.StreamHandler(),  # Log to the console
            logging.FileHandler('multiplanetaryListUpdBot.log', mode='a')  # Log to a file (append mode)
            ]
        )
    warnings.filterwarnings('ignore', category=AstropyWarning)
    sqliteCursor.execute("DROP TABLE IF EXISTS stars")
    with sqliteConn:
        sqliteConn.execute("""
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
          planets INTEGER,
          altNames TEXT
        );
      """)

    sqliteCursor.execute("DROP TABLE IF EXISTS simbad")
    with sqliteConn:
        sqliteConn.execute("""
         CREATE TABLE simbad (
         name TEXT,
         ra TEXT,
         dec TEXT,
         mag REAL,
         dist REAL,
         ids TEXT,
         PRIMARY KEY (ra, dec)
        );
      """)
   
    print("Retrieving data from Exoplanet ...")
    getDataFromExoplanet(None)    
    #getDataFromExoplanet("exoplanet.csv")
    print("Retrieving data from Simbad ...")
    getDataFromSimbadSite()
    #print("Retrieving data from NASA ...")
    #getDataFromNASA(None)                     
    #getDataFromNASA("NASA.out")
    print("Retrieving data from Wikipedia ...")
    getDataFromWikipedia(None)
    #getDataFromWikipedia("wiki.out")
    print("Generating 'tabella.wiki' ...")
    generateWikitable("tabella.wiki")
    print("Done.")
    print("Copy the content of the file 'tabella.wiki' in the correct place inside https://it.wikipedia.org/wiki/Sistemi_multiplanetari. Check the result before publishing!")


if __name__ == "__main__":
    main()
