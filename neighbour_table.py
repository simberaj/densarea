import common
import arcpy

TOLERANCE = '1 Centimeters'
EXTERIOR_ID = -1

def table(zones, idFld, output, exterior=True, selfrel=True):
  common.debug('running neighbour table', zones, idFld, output, exterior, selfrel)
  with common.PathManager(output) as pathman:
    if exterior:
      common.progress('mapping zone surroundings')
      buffer = pathman.tmpFC()
      arcpy.Buffer_analysis(zones, buffer, '50 Meters', dissolve_option='ALL')
      common.progress('creating exterior zone')
      erased = pathman.tmpFC()
      arcpy.Erase_analysis(buffer, zones, erased, TOLERANCE)
      common.progress('identifying exterior zone')
      common.calcField(erased, idFld, EXTERIOR_ID,
        common.pyTypeOfField(zones, idFld))
      # common.progress('eliminating sliver polygons')
      common.progress('merging exterior zone')
      jointo = pathman.tmpFC()
      arcpy.Merge_management([zones, erased], jointo)
    else:
      jointo = zones
    common.progress('finding neighbours')
    joined = pathman.tmpFC()
    fm = arcpy.FieldMappings()
    fm.addFieldMap(common.fieldMap(zones, idFld, common.NEIGH_FROM_FLD, 'FIRST'))
    fm.addFieldMap(common.fieldMap(jointo, idFld, common.NEIGH_TO_FLD, 'FIRST'))
    arcpy.SpatialJoin_analysis(zones, jointo, joined, 'JOIN_ONE_TO_MANY', 'KEEP_COMMON', fm, 'INTERSECT', TOLERANCE)
    common.progress('converting to neighbour table')
    fm2 = arcpy.FieldMappings()
    fm.addFieldMap(common.fieldMap(joined, common.NEIGH_FROM_FLD, common.NEIGH_FROM_FLD, 'FIRST'))
    fm.addFieldMap(common.fieldMap(joined, common.NEIGH_TO_FLD, common.NEIGH_TO_FLD, 'FIRST'))
    if selfrel:
      query = ''
    else:
      query = common.safeQuery('[{}] <> [{}]'.format(
        common.NEIGH_FROM_FLD, common.NEIGH_TO_FLD), joined)
    arcpy.TableToTable_conversion(joined, pathman.getLocation(), pathman.getOutputName(), query, fm2)
    common.clearFields(output, [common.NEIGH_FROM_FLD, common.NEIGH_TO_FLD])
  return output

if __name__ == '__main__':
  with common.runtool(5) as parameters:
    zones, idFld, output, exteriorStr, selfrelStr = parameters
    exterior = common.toBool(exteriorStr, 'exterior relationship record switch')
    selfrel = common.toBool(selfrelStr, 'self-neighbourhood record switch')
    table(zones, idFld, output, exterior, selfrel)

