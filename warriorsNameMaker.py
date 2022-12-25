#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import random
prefixes =['forest', 'moor', 'oak', 'bracken', 'elk', 'sheep', 'flame', 'night', 'pine', 'fish']
suffixes = ['claw', 'paw', 'paw', 'kit', 'kit', 'star', 'kit', 'star', 'paw', 'kit']
colors = ['orange', 'black', 'white']
causeofdeath =['enemycats','rouge','illness']
numberofkits = [1,3,0]
eyecolor = ['blue','yellow','green']
patterns = ['stipes','bloches','plain']
clan = ['thunder','wind','river']
fluffyness = ['notvery','OMGFLUFFY','kindafluff']
patterncolor = ['orange black','grey','yellow']
gender = ['she-cat','tom']
nuberofapps = [1,2,0]



prefixes.append('fox')
suffixes.append('tail')
prefixes.append('squirrel')
prefixes.append('bird')
prefixes.append('sparrow')
prefixes.append('chimpmunk')
prefixes.append('rabbit')
prefixes.append('vole')
prefixes.append('hare')
prefixes.append('olive')
prefixes.append('tall')
prefixes.append('leaf')
prefixes.append('arrow')
prefixes.append('daisy')
prefixes.append('yarrow')
prefixes.append('yellow')
prefixes.append('claw')
prefixes.append('swamp')
prefixes.append('scar')
prefixes.append('long')
prefixes.append('mose')
prefixes.append('mole')
prefixes.append('horse')
prefixes.append('brick')
suffixes.append('stripe')
suffixes.append('berrie')
suffixes.append('eye')
suffixes.append('pool')
suffixes.append('moon')
suffixes.append('ear')
suffixes.append('twig')
suffixes.append('glass')
suffixes.append('nose')
suffixes.append('muzzle')
suffixes.append('prey')
suffixes.append('strike')
suffixes.append('hunt')
suffixes.append('fur')
suffixes.append('fur')
suffixes.append('nut')
suffixes.append('pelt')
suffixes.append('flower')
suffixes.append('thorn')
suffixes.append('bush')
suffixes.append('leaf')
suffixes.append('whisker')
suffixes.append('heart')
suffixes.append('blaze')
suffixes.append('bone')
prefixes.append('dog')
prefixes.append('grow')
prefixes.append('vine')
prefixes.append('mint')
prefixes.append('leopard')
prefixes.append('nettle')
prefixes.append('mole')
prefixes.append('orange')
prefixes.append('badger')
prefixes.append('sky')
prefixes.append('noble')
prefixes.append('jay')
prefixes.append('starling')
prefixes.append('bedraggled')
prefixes.append('ember')
prefixes.append('sun')
prefixes.append('silver')
prefixes.append('storm')
prefixes.append('gold')
prefixes.append('bronze')
prefixes.append('lake')
prefixes.append('moth')
prefixes.append('warm')
prefixes.append('meadow')
prefixes.append('rat')
prefixes.append('glow')
prefixes.append('seal')
prefixes.append('egg')
prefixes.append('bright')
prefixes.append('bare')
prefixes.append('run')
prefixes.append('cold')
prefixes.append('hot')
prefixes.append('chill')
prefixes.append('coal')
prefixes.append('ash')
prefixes.append('dot')
prefixes.append('spotted')
prefixes.append('broke')
prefixes.append('feather')
prefixes.append('jagged')
prefixes.append('tree')
prefixes.append('moon')
prefixes.append('orange')
prefixes.append('light')
prefixes.append('shimmer')
prefixes.append('spot')
prefixes.append('gleam')
causeofdeath.append('badger')
causeofdeath.append('fox')
causeofdeath.append('burning')
causeofdeath.append('falling')
causeofdeath.append('drowning')
causeofdeath.append('sufficating')
causeofdeath.append('choking')
causeofdeath.append('poisoned')
eyecolor.append('blind')
eyecolor.append('light green')
eyecolor.append('dark blue')
eyecolor.append('pale yellow')
eyecolor.append('light blue')
eyecolor.append('amber')
eyecolor.append('light brown')
eyecolor.append('black')
eyecolor.append('grey')
eyecolor.append('dark green')
eyecolor.append('light orange')
eyecolor.append('blue and green')
eyecolor.append('light green and dark blue')
eyecolor.append('brown and black')


numberOfPrefixes = len(prefixes)
numberOfSuffixes = len(suffixes)
totalNumberOfPossibleNames = numberOfPrefixes*numberOfSuffixes



randomPrefix = random.choice(prefixes)
randomSuffix = random.choice(suffixes)
#print(numberOfPrefixes,' Number of Prefixes ',numberOfSuffixes, ' Number Of Suffixes ',totalNumberOfPossibleNames,' Total # Possible Names')
print('Your name is ', randomPrefix+randomSuffix, ' you are in ',random.choice(clan) ,' clan')
print('You are a ', random.choice(colors), ' ', random.choice(gender), ' cat with ', random.choice(patterncolor),' ',random.choice(patterns),'and your eye color is ',random.choice(eyecolor) )

print('your fluffyness level is ', random.choice(fluffyness) )


if not (randomSuffix == 'kit' or randomSuffix == 'paw'):
    totalNumberOfApprenntices = random.choice(nuberofapps)

    if totalNumberOfApprenntices > 0:
        print("You have ",totalNumberOfApprenntices,' apprentices, their names are: ')
        for thisApprentice in range(0,totalNumberOfApprenntices):
            apprenticeName = random.choice(prefixes)+random.choice(suffixes)
            print(apprenticeName)

    totalnumberofkits = random.choice(numberofkits)

    if totalnumberofkits > 0:
        print("You have ",totalnumberofkits,' kits, their names are: ')
        for thiskit in range(0,totalnumberofkits):
            kitsName = random.choice(prefixes)+random.choice(suffixes)
            print(kitsName)




print('you died of ', random.choice(causeofdeath) )


def getName():
    randomPrefix = random.choice(prefixes)
    randomSuffix = random.choice(suffixes)
    fullName = randomPrefix+randomSuffix
    return fullName

# EXAMPLES
#Your name is  skyblaze  you are in  wind  clan
#You are a  white   she-cat  cat with  orange black   bloches and your eye color is  blue
#You have  2  apprentices, their names are:
#pinekit
#treebush
#You have  3  kits, their names are:
#fishbush
#eggfur
#runmoon
#you died of  falling
#
#
#
#Your name is  clawkit  you are in  thunder  clan
#You are a  black   tom  cat with  orange black   stipes and your eye color is  blue
#your fluffyness level is  OMGFLUFFY
#you died of  rouge
#
#
#
#Your name is  molepelt  you are in  wind  clan
#You are a  black   tom  cat with  orange black   stipes and your eye color is  blue and green
#your fluffyness level is  OMGFLUFFY
#You have  1  apprentices, their names are:
#nettleeye
#You have  3  kits, their names are:
#volestripe
#starlingear
#leopardkit
#you died of  fox
