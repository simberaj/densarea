import common
import loaders
import regional
import operator

DENSITY_COEF = 1000000
DENS_CALLER = operator.itemgetter('dens')

class DensityZone(regional.Zone):
  def __init__(self, id, **kwargs):
    regional.Zone.__init__(self, id, **kwargs)
    self._set('dens', self.get('mass') / float(self.get('area')) * DENSITY_COEF)
  
  
class DensityAreal(regional.Region):
  cached = ['mass', 'area']

  def __init__(self, zone):
    regional.Region.__init__(self, zone.getID())
    self.bind(zone)
  
  def relabel(self):
    self.setID(max(self.zones, key=regional.MASS_CALLER).getID())
    
  def _modified(self):
    self._set('dens', (self.get('mass') / float(self.get('area')) * DENSITY_COEF) if self else 0.0)
    
  def isAccepted(self, zone, thrDens):
    return (self.get('mass') + zone.get('mass')) / (self.get('area') + zone.get('area')) * DENSITY_COEF >= thrDens
  
  def getNextZone(self, doMergeEnclaves=False):
    neighs = self.getNeighZones()
    while neighs:
      if len(neighs) == 1: return neighs.pop()
      densest = max(neighs, key=DENS_CALLER)
      if doMergeEnclaves:
        enclzones = self.potentialEnclaves([densest])
        if enclzones:
          densWithEncl = (sum(zone.get('mass') for zone in enclzones) / 
              sum(zone.get('area') for zone in enclzones) * DENSITY_COEF)
          neighs.remove(densest)
          if densWithEncl >= max(neighs, key=DENS_CALLER):
            return densest
        else:
          return densest
      else:
        return densest
    return None

  
          
def createAreals(zones, thrDens):
  areals = []
  common.progress('creating areals')
  for zone in zones.values():
    if zone.get('dens') >= thrDens:
      areals.append(DensityAreal(zone))
  return areals

def regionalize(zones, thrDens, minPop, doMergeEnclaves=True):
  areals = createAreals(zones, thrDens)
  todo = set(areals)
  while todo:
    areal = todo.pop()
    # if areal.getID().startswith('554782'): common.debug('working on', areal)
    while True:
      nextZone = areal.getNextZone(doMergeEnclaves)
      # if areal.getID().startswith('554782'): common.debug('next zone for', areal, nextZone)
      if nextZone is None:
        break
      elif nextZone.isAssigned(): # two areals connected, merge them
        other = nextZone.getRegion()
        todo.discard(other)
        areal.merge(other)
      elif areal.isAccepted(nextZone, thrDens):
        # if areal.getID().startswith('554782'): common.debug('adding', areal, nextZone)
        areal.bind(nextZone)
        if doMergeEnclaves:
          areal.includeEnclaves()
      else:
        break
      # if areal.getID().startswith('554782'): common.debug('state', areal, areal.get('dens'), areal.get('mass'))
  mergeAdjacent(areals)
  # if doMergeEnclaves:
    # areals = resolveEnclaves(areals, thrDens, minPop)
  eraseSmall(areals, minPop)
  for areal in areals:
    if areal:
      areal.relabel()
  # if doMergeEnclaves:
    # warnUnderdens(areals, thrDens)
    
def warnUnderdens(areals, thrDens):
  for areal in areals:
    if areal and areal.get('dens') < thrDens:
      common.warning('could not densify areal {}, leaving density at {}'.format(areal.getID(), areal.get('dens')))

def eraseSmall(areals, minPop):
  for areal in areals:
    # common.debug(areal, areal.get('mass'), minPop)
    if areal.get('mass') < minPop:
      areal.erase()
  
def mergeAdjacent(areals):
  for areal in areals:
    if areal:
      for cont in areal.getNeighRegions():
        areal.merge(cont)
          
  
def resolveEnclaves(areals, thrDens, minMass=0):
  common.progress('resolving enclaves')
  newAreals = []
  while areals:
    areal = areals.pop()
    if areal:
      if areal.isInEnclave():
        areal.erase()
      else:
        areal.includeEnclaves() # include all pockets
        newAreals.append(areal)
        areals.extend(areal.densify(thrDens, minMass))
  return newAreals

def groupDensity(zonelist):
  mass = 0.0
  area = 0.0
  for zone in zonelist:
    mass += zone.get('mass')
    area += zone.get('area')
  return mass / area
      
def delimitDensityAreals(zones, idFld, popFld, thrDens, minPop, targetFld, neighTable=None, doMergeEnclaves=True):
  common.progress('loading areal data')
  loader = loaders.RegionalLoader()
  # common.progress('calculating zone densities')
  areaFld = common.ensureShapeAreaField(zones)
  inSlots = {'id' : idFld, 'mass' : popFld, 'area' : areaFld}
  loader.sourceOfZones(zones, inSlots, targetClass=DensityZone)
  loader.possibleNeighbourhood(neighTable, exterior=True)
  loader.load()
  zones = loader.getZoneDict()
  common.progress('delimiting areals')
  regionalize(zones, thrDens, minPop, doMergeEnclaves)
  common.progress('saving data')
  loader.addZoneOutputSlot('assign', targetFld, require=True)
  loader.outputZones(zones)

if __name__ == '__main__':
  with common.runtool(8) as parameters:
    parameters[3] = common.toFloat(parameters[3], 'threshold density')
    parameters[4] = common.toFloat(parameters[4], 'minimum areal population')
    parameters[7] = common.toBool(parameters[7], 'enclave merge switch')
    # import cProfile
    # cProfile.run('delimitDensityAreals(*parameters)')
    delimitDensityAreals(*parameters)