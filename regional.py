from collections import defaultdict
import operator
import itertools
# import common # only for debug

MASS_CALLER = operator.itemgetter('mass')

class RegionalUnit:
  '''A measurable regional unit - a pseudoabstract superclass of Zone and Region allowing both of them to provide IDs.'''
  id = None

  def __init__(self, id, **kwargs):
    self.id = id
    self._props = kwargs

  def _set(self, key, value):
    self._props[key] = value
  
  def get(self, key):
    return self._props[key]
  
  def __getitem__(self, key):
    return self._props[key]
  
  def setID(self, id):
    self.id = id
    
  def getID(self):
    return self.id

  # for interaction calculations (excluding zones by their region)
  def getRegion(self):
    return None
    
  def __repr__(self):
    return '<{} {}{}>'.format(self.__class__.__name__, self.id,
        ' (' + '|'.join(unicode(self.get(key)) for key in self.display) + ')' if self.display else '')
  
    
class Neighbour:
  '''A simple superclass allowing neighbourhood formalization.'''
  def __init__(self):
    self.neighbours = set()
  
  def addNeighbour(self, zone):
    self.neighbours.add(zone)
  
  def hasNeighbour(self, zone):
    return (zone in self.neighbours)
  
  def setNeighbours(self, neighs):
    if neighs:
      self.neighbours.update(neighs)
  
  def getNeighbours(self):
    return self.neighbours

class Exterior(Neighbour):
  def getID(self):
    return -1

  def __repr__(self):
    return '<Exterior>'

exterior = Exterior()
    
class Zone(RegionalUnit, Neighbour):
  delegation = 'region'
  coreable = True
  display = []
  
  def __init__(self, id, **kwargs):
    RegionalUnit.__init__(self, id, **kwargs)
    Neighbour.__init__(self)
    self.clearRegion()
  
  def setRegion(self, region):
    self.region = region
  
  def getRegion(self):
    return self.region
  
  def clearRegion(self):
    self.region = None
  
  def getRegionID(self):
    if self.region is not None:
      return self.region.getID()
    else:
      return None
  
  def isAssigned(self):
    return self.region is not None
  
  def isInRegion(self, region):
    return self.region is region

  def getNeighRegions(self):
    contig = set(zone.getRegion() for zone in self.neighbours if zone is not exterior)
    contig.discard(None)
    return contig
  
  def hasNeighRegion(self, region):
    for zone in self.neighbours:
      if zone is not exterior and region is zone.getRegion():
        return True
    return False
  
  def isOnRegionEdge(self):
    # zone will be on edge if any neighbour's region is different from self
    for neigh in self.neighbours:
      if neigh is exterior or not neigh.isInRegion(self.region):
        return True
    else:
      return False



class Region(RegionalUnit):
  cached = []
  display = ['count']

  def __init__(self, id):
    RegionalUnit.__init__(self, id)
    self._cache = self._props
    self.clear()
    
  def __nonzero__(self):
    '''The region is True if it contains any zone (any assignment).'''
    return bool(self.zones)
  
  def add(self, zone):
    self.zones.add(zone)
    self.cuts = None
    self._addToCache(zone)
  
  def remove(self, zone):
    self.zones.remove(zone)
    self.cuts = None
    self._subFromCache(zone)
    
  def _addToCache(self, zone):
    for key in self.cached:
      self._cache[key] += zone.get(key)
    self._cache['count'] += 1
    self._modified()
  
  def _subFromCache(self, zone):
    for key in self.cached:
      self._cache[key] -= zone.get(key)
    self._cache['count'] -= 1
    self._modified()
  
  def _modified(self):
    pass
    
  def clear(self):
    self.zones = set()
    self.cuts = None
    for key in self.cached:
      self._cache[key] = 0
    self._cache['count'] = 0
    
  def bind(self, zone):
    zone.setRegion(self)
    self.add(zone)

  def unbind(self, zone):
    zone.clearRegion()
    self.remove(zone)
  
  def detach(self):
    for zone in self.zones:
      zone.clearRegion()
  
  def erase(self):
    zonelist = list(self.zones)
    while zonelist:
      self.unbind(zonelist.pop())
  
  def merge(self, reg):
    otherZones = list(reg.getZones())
    reg.erase()
    for zone in otherZones:
      self.bind(zone)
  
  def getZones(self):
    return self.zones
  
  def getCuts(self):
    if self.cuts is None:
      self.cuts = self._calcCuts()
    return self.cuts
  
  def _calcCuts(self):
    '''Calculates region cut points - zones that would cause some other zones of the region to become exclaves.'''
    if not self.zones: return {}
    artic = defaultdict(list) # articulation points and portions they hide from the start zone
    togo = [] # DFS stack
    ins = {} # enter time of DFS
    lows = {} # lowpoint function
    neighs = {} # neighbourhood list
    children = defaultdict(list) # DFS tree children of given zones
    # find start zone and load neighbour matrix
    for zone in self.zones:
      if not togo: togo.append(zone)
      neighs[zone] = [neigh for neigh in zone.getNeighbours() if
          neigh is not exterior and neigh.isInRegion(self)]
    # start DFS
    root = togo[-1]
    counter = 0
    while togo:
      now = togo[-1]
      if now not in ins: # enter vertex
        ins[now] = counter
        lows[now] = counter
        counter += 1
      for neigh in neighs[now]: # add all unvisited neighbours to stack
        if neigh not in ins:
          children[now].append(neigh)
          togo.append(neigh)
          break
      else: # exit vertex (if no unvisited neighbours found)
        if now is not root: # root has different articulation procedures
          counter += 1
          for neigh in neighs[now]: # all have been visited; update lowpoint
            if lows[neigh] < lows[now]:
              lows[now] = lows[neigh]
          for neigh in neighs[now]: # check if neigh is articulation
            if lows[neigh] >= ins[now]:
              for sub in artic[now]:
                if neigh in sub:
                  break
              else:
                artic[now].append(self.subtree(children, neigh)) # get what neigh separates from now
            elif now not in children[neigh] and lows[now] > ins[neigh]:
              lows[now] = ins[neigh]
        del togo[-1]
    for child in children[root][1:]: # if root has 2+ children, it is articulation
      artic[now].append(self.subtree(children, child))
    # if artic:
      # common.debug('articulation of %s detected: %s' % (self, artic))
    return artic
  
  @staticmethod
  def subtree(children, root):
    stack = [root]
    tree = set()
    while stack:
      now = stack.pop()
      tree.add(now)
      stack.extend(children[now])
    return tree
  
  def getNeighZones(self, includeExterior=False):
    contig = set()
    for zone in self.zones:
      contig.update(zone.getNeighbours())
    contig.difference_update(self.zones)
    if not includeExterior:
      contig.discard(exterior)
    return list(contig)
  
  def getNeighRegions(self):
    return list(set(zone.getRegion() for zone in self.getNeighZones()) - set([None]))
    
  def includeEnclaves(self):
    for enclave in self.enclaves():
      # common.debug(self, 'enclave found', enclave)
      for zone in enclave:
        self.bind(zone)
  
  def isInEnclave(self):
    '''Returns True if the region is entirely enclosed with an another
    region, allowing for unassigned space in between.'''
    todo = set(self.getNeighZones(includeExterior=True))
    if exterior in todo:
      # common.debug(self, 'exterior contact: direct')
      return False
    visited = set(self.zones)
    found = set([self])
    while todo:
      current = todo.pop()
      if current is exterior:
        # common.debug(self, 'exterior contact: indirect')
        return False # never an enclave when exterior touched
      elif current.isAssigned():
        found.add(current.getRegion())
        if len(found) > 2: # at least 2 other regions found
          # common.debug(self, 'diversity contact:', found)
          return False
      else:
        todo.update(current.getNeighbours().difference(visited))
      visited.add(current)
    # common.debug(self, 'in an enclave')
    return True
  
  def enclaves(self, additional=[]):
    startPoints = set(zone for zone in self.getNeighZones() if not zone.isAssigned())
    return self._enclavesearch(startPoints, additional)
  
  def _enclavesearch(self, start, addblock=[]):
    free = set()
    block = set(self.zones)
    block.update(addblock)
    encl = set()
    while start:
      zone = start.pop()
      found, tree = self.searchTree(zone, block)
      if found:
        free.update(tree)
        start.difference_update(tree)
      else:
        encl.add(tuple(tree))
    return list(encl)
    
  def potentialEnclaves(self, additional):
    return list(set(itertools.chain.from_iterable(self._enclavesearch(self._potentialEnclaveCandidates(additional), additional))))
  
  def _potentialEnclaveCandidates(self, additional):
    stp = set()
    for zone in additional:
      stp.update(neigh for neigh in zone.getNeighbours() if neigh is not exterior and not neigh.isAssigned())
    return stp
        
  def enclaveZones(self, additional=[]):
    return list(set(itertools.chain.from_iterable(self.enclaves(additional))))
      
  def searchTree(self, start, block):
    found = False
    stack = [start]
    tree = set(stack)
    while stack and not found:
      current = stack.pop()
      for neigh in current.getNeighbours():
        if neigh not in block and neigh not in tree:
          if neigh is exterior or neigh.isAssigned():
            found = True # region cannot be self because self are in block
            break
          else: # not assigned, see through it
            stack.append(neigh)
            tree.add(neigh)
    return found, tree
    
  def connectedComponents(self):
    notfound = set(self.zones)
    visited = set()
    comps = []
    while notfound:
      stack = set([notfound.pop()])
      comps.append([])
      while stack:
        current = stack.pop()
        comps[-1].append(current)
        visited.add(current)
        notfound.discard(current)
        stack.update(zone for zone in current.getNeighbours() if 
            zone is not exterior and zone.isInRegion(self) and zone not in visited)
    assert sum(len(comp) for comp in comps) == len(self.zones)
    return comps