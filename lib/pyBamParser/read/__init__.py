#Dan Blankenberg
import struct
# import logging

from ..util.packer import pack_int8, unpack_int8, pack_uint8, unpack_uint8, pack_int16, unpack_int16, pack_uint16, unpack_uint16, pack_int32, unpack_int32, pack_uint32, unpack_uint32, pack_int64, unpack_int64, pack_uint64, unpack_uint64
from ..util import NULL_CHAR

BAM_NO_QUAL = 0xFF #255
SEQ_4_BIT_TO_SEQ = list( '=ACMGRSVTWYHKDBN' )
# OP_INTS_TO_CHARS = {0:'M', 1:'I', 2:'D', 3:'N', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X'}

CIGAR_OP = list( 'MIDNSHP=X' )

NULL_TERMINATED_TAGS = [ 'Z', 'H' ]
TAG_TYPE_TO_VALUE_LENGTH = {
                            'A': 1,
                            'c': 1,
                            'C': 1,
                            's': 2,
                            'S': 2,
                            'i': 4,
                            'I': 4,
                            'f': 4,
                            'Z': 1,
                            'H': 2,
                            'B': 1
                            }
TAG_TYPE_TO_STRUCT_TYPE = {
                            'A': 'c',
                            'c': 'b',
                            'C': 'B',
                            's': 'h',
                            'S': 'H',
                            'i': 'i',
                            'I': 'I',
                            'f': 'f',
                            'Z': 'c',
                            'H': 'c',
                            'B': 'c'
                            }

TAG_TYPE_TO_SAM_TYPE = {
                            'A': 'A',
                            'c': 'i',
                            'C': 'i',
                            's': 'i',
                            'S': 'i',
                            'i': 'i',
                            'I': 'i',
                            'f': 'f',
                            'Z': 'Z',
                            'H': 'H',
                            'B': 'B'
                            }

BAM_READ_BEGIN_UNPACKER = struct.Struct( "<iiIIiiii" ).unpack
READ_GROUP_RECORD_TAG = 'RG'

SEQ_UNPACKERS = {} #cache seq unpackers for reuse
QUAL_UNPACKERS = {}
CIGAR_UNPACKERS = {}

#cache seq bit to 4 bit lookup
SEQ_4_BIT_TO_SEQ_4 = [ SEQ_4_BIT_TO_SEQ[ i >> 4 ]  for i in range( 256 ) ]
SEQ_4_BIT_TO_SEQ_F = [ SEQ_4_BIT_TO_SEQ[ i & 0x0f ]  for i in range( 256 ) ]



#TODO: FIXME: make parsing occur on-demand, by attribute, not at initialization
class BAMRead( object ):
    
    def __init__( self, data, reader ): #TODO: instead of reader, can we take list of sequences? can we omit block size?
        self._block_size = len( data )
        self._reader = reader
        self._ref_id, self._pos, self._bin_mq_nl, self._flag_nc, self._l_seq, self._next_ref_id, self._next_pos, self._t_len = BAM_READ_BEGIN_UNPACKER( data[ :32 ] )
        self._bin = self._bin_mq_nl >> 16
        self._mapq = self._bin_mq_nl >> 8 & 0xff
        self._l_read_name = self._bin_mq_nl & 0xff
        self._flag = self._flag_nc >> 16
        self._n_cigar_op = self._flag_nc & 0xff

        self.__data = data
        self.__not_parsed_1 = True
        self.__not_parsed_2 = True
        self.__not_parsed_3 = True
        self.__not_parsed_4 = True
        self.__not_parsed_5 = True
        
        self.__zero_based_end_position = None
        self.__seq_string = None
        self.__qual_list = None
        self.__is_seq_reverse_complement = None
        self.__read_group = None
        self.__read_group_parsed = False
        self.__reference_name = None
        self.__reference = None
        self.__contiguous_blocks = None
        self.__indels = None

    def __parse_block_1( self ):
        if self.__not_parsed_1:
            self.__not_parsed_1 = False
            self._block_offset = 32 + self._l_read_name
            self._read_name = self.__data[ 32:self._block_offset].rstrip( NULL_CHAR )
            
    def __parse_block_2( self ):
        if self.__not_parsed_2:
            self.__parse_block_1()
            self.__not_parsed_2 = False
            n_cigar_op = self._n_cigar_op
            cigar_op_len = 4 * n_cigar_op
            self._new_block_offset = self._block_offset + cigar_op_len
            cigar_unpacker = CIGAR_UNPACKERS.get( n_cigar_op, None )
            if cigar_unpacker is None:
                cigar_unpacker = struct.Struct( "<" +"I" * n_cigar_op ).unpack
                CIGAR_UNPACKERS[ n_cigar_op ] = cigar_unpacker
            self._cigar = cigar_unpacker( self.__data[ self._block_offset:self._new_block_offset ] )
            self._cigar_list = []
            for cigar in self._cigar:
                op_len = cigar >> 4
                op = cigar & 0x07
                self._cigar_list.append( ( op_len, op ) )
            
    def __parse_block_3( self ):
        if self.__not_parsed_3:
            self.__parse_block_2()
            self.__not_parsed_3 = False
            seq_bin_len = ( self._l_seq + 1 ) / 2
            self._block_offset = self._new_block_offset
            self._new_block_offset += seq_bin_len
            seq_unpacker = SEQ_UNPACKERS.get( seq_bin_len, None )
            if seq_unpacker is None:
                seq_unpacker = struct.Struct( "<" + "B" * seq_bin_len ).unpack
                SEQ_UNPACKERS[ seq_bin_len ] = seq_unpacker #cache this unpacker for use later
            self._seq = seq_unpacker( self.__data[ self._block_offset: self._new_block_offset ] )
            half_seq_len = self._l_seq / 2
            seq = []
            for i, s in enumerate( self._seq, start=1 ):
                seq.append( SEQ_4_BIT_TO_SEQ_4[ s ] )
                if i <= half_seq_len:
                    seq.append( SEQ_4_BIT_TO_SEQ_F[ s ] )
                else:
                    break
            self._seq_list = seq
            self._block_offset = self._new_block_offset
            self._new_block_offset += self._l_seq
            
    def __parse_block_4( self ):
        if self.__not_parsed_4:
            self.__parse_block_3()
            self.__not_parsed_4 = False
            qual_unpacker = QUAL_UNPACKERS.get( self._l_seq, None )
            if qual_unpacker is None:
                qual_unpacker = struct.Struct( "<" + "c" * self._l_seq ).unpack
                QUAL_UNPACKERS[ self._l_seq ] = qual_unpacker
            self._qual = map( ord, qual_unpacker( self.__data[ self._block_offset: self._new_block_offset ] ) )
            if self._qual[0] == BAM_NO_QUAL:
                self._qual = None
            
    def __parse_block_5( self ):
        if self.__not_parsed_5:
            self.__parse_block_4()
            self.__not_parsed_5 = False
            self._aux_data = []
            while self._block_size > self._new_block_offset:
                self._block_offset = self._new_block_offset
                self._new_block_offset += 2
                tag = self.__data[ self._block_offset: self._new_block_offset ]
                self._block_offset = self._new_block_offset
                self._new_block_offset += 1
                val_type = self.__data[ self._block_offset: self._new_block_offset ]
                if val_type in NULL_TERMINATED_TAGS:
                    value = []
                    while True:
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 1
                        val = self.__data[ self._block_offset: self._new_block_offset ]
                        if val == NULL_CHAR:
                            break
                        else:
                            value.append( val )
                    value = ( "".join( value ), )
                else:
                    if val_type == 'B':
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 1
                        val_type = self.__data[ self._block_offset: self._new_block_offset ]
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 4
                        tag_length = unpack_int32( self.__data[ self._block_offset: self._new_block_offset ] )[ 0 ]
                    else:
                        tag_length = 1
                    val_size = TAG_TYPE_TO_VALUE_LENGTH[ val_type ] * tag_length
                    self._block_offset = self._new_block_offset
                    self._new_block_offset += val_size
                    value = self.__data[ self._block_offset: self._new_block_offset ]
                    value = struct.unpack( "<" + TAG_TYPE_TO_STRUCT_TYPE[ val_type ] * tag_length, value )
                self._aux_data.append( ( tag, val_type, value ) )
            self.__data = None
    
    def get_read_name( self ):
        self.__parse_block_1()
        return self._read_name
    def _get_bam_read_name( self ):
        return self.get_read_name() + NULL_CHAR
    
    def get_flag( self ):
        return self._flag
    
    def get_reference( self ):
        if self.__reference is None:
            self.__reference = self._reader.get_reference_by_id( self._ref_id )
        return self.__reference
        
    def get_reference_name( self  ):
        if self.__reference_name is None:
            self.__reference_name = self._reader.get_reference_name_by_id( self._ref_id )
        return self.__reference_name
    
    def get_reference_id( self ):
        return self._ref_id
    
    def _get_bam_ref_id( self ):
        return pack_int32( self._ref_id )
    
    def get_rnext( self ):
        return self._reader.get_reference_by_id( self._next_ref_id )
    
    def get_rnext_name( self ):
        return self._reader.get_reference_name_by_id( self._next_ref_id, self._ref_id  )
    
    def _get_bam_rnext_id( self ):
        return pack_int32( self._next_ref_id )
    
    def get_position( self, one_based=True ):
        return self._pos + one_based
    def get_position_zero_based( self ):
        return self._pos
    def _get_bam_pos( self ):
        return pack_int32( self.get_position( False ) )

    def to_ref_coord( self, read_coord, one_based=True ):
        blocks = self.get_contiguous_blocks( one_based=one_based )
        return self._to_ref_coord( blocks, read_coord )

    @staticmethod
    def _to_ref_coord( blocks, read_pos ):
        for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
            if direction == 1:
                hit = read_start <= read_pos < read_end
            elif direction == -1:
                hit = read_end < read_pos <= read_start
            if hit:
                return direction * read_pos + offset
        # logging.warn('No hit on read coordinate {}.'.format(read_pos))

    #TODO: def _to_read_coord( self, blocks, ref_pos ):

    def get_contiguous_blocks( self, one_based=True ):
        """Return a list of blocks of aligned bases.
        Each block is a 6-tuple:
        1. The start of the block, in read coordinates.
        2. The end of the block, in read coordinates.
        3. The start of the block, in reference coordinates.
        4. The end of the block, in reference coordinates.
        5. The offset between read and reference coordinates (offset = ref - direction * read).
        6. The direction of the read (1 for forward, -1 for reverse).
        """
        # Generate contiguous blocks, if none is cached.
        if self.__contiguous_blocks is None:
            cigar = self.get_cigar()
            ref_pos = self.get_position( one_based=one_based )
            reverse = self.is_seq_reverse_complement()
            if reverse:
                read_len = self.get_l_seq()
            else:
                read_len = None
            self.__contiguous_blocks = self._get_contiguous_blocks( ref_pos, cigar, reverse, read_len )
        return self.__contiguous_blocks

    @staticmethod
    def _get_contiguous_blocks( ref_pos, cigar, reverse, read_len ):
        """Do the actual CIGAR string parsing to generate the list of aligned bases."""
        ref_pos_start = ref_pos
        if reverse:
            direction = -1
            read_pos = read_len
        else:
            direction = 1
            read_pos = 1
        read_pos_start = read_pos
        blocks = []
        # logging.info('Ref starting at {}, read {}.'.format(ref_pos, read_pos))
        while cigar:
            cigar_size, cigar_op = cigar.pop( 0 )
            # logging.info('Saw {}{}.'.format(cigar_size, OP_INTS_TO_CHARS[cigar_op]))
            # 0: M alignment match (can be a sequence match or mismatch)
            # 7: = sequence match
            # 8: X sequence mismatch
            if cigar_op in ( 0, 7, 8 ):
                ref_pos += cigar_size
                read_pos += cigar_size * direction
            # 1: I insertion
            # 4: S soft clipping (clipped sequences present in SEQ)
            elif cigar_op in ( 1, 4 ):
                offset = ref_pos - direction * read_pos
                if read_pos_start != read_pos:
                    blocks.append( ( read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction ) )
                read_pos += cigar_size * direction
                read_pos_start = read_pos
                ref_pos_start = ref_pos
            # 2: D deletion from the reference
            # 3: N skipped region from the reference
            elif cigar_op in ( 2, 3 ):
                offset = ref_pos - direction * read_pos
                blocks.append( ( read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction ) )
                ref_pos += cigar_size
                read_pos_start = read_pos
                ref_pos_start = ref_pos
            # 5: H hard clipping (clipped sequences NOT present in SEQ)
            # 6: P padding (silent deletion from padded reference)
            elif cigar_op in ( 5, 6 ):
                pass
            else:
                pass #logging.warn('unknown cigar_op {} {}'.format(cigar_op, cigar_size))
            # logging.info('Ref now {}, read {}.'.format(ref_pos, read_pos))
        offset = ref_pos - direction * read_pos
        if read_pos_start != read_pos:
            blocks.append( ( read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction ) )
        return blocks

    def indel_at( self, position, check_insertions=True, check_deletions=True, one_based=True ):
        """Does the read contain an indel at the given position?
        Return True if the read contains an insertion at the given position
        (position must be the base before the insertion event) or if the read
        contains a deletion where the base at position is deleted. Return False
        otherwise."""
        insertions, deletions = self.get_indels( one_based=one_based )
        return self._indel_at( position, insertions, deletions, check_insertions=check_insertions,
                               check_deletions=check_deletions )

    @staticmethod
    def _indel_at( position, insertions, deletions, check_insertions=True, check_deletions=True ):
        if check_insertions:
            for insertion in insertions:
                if insertion[0] == position:
                    return True
        if check_deletions:
            for deletion in deletions:
                if deletion[0] < position < deletion[0] + deletion[1] + 1:
                    return True
        return False

    def get_indels( self, one_based=True ):
        """Return a data structure containing all indels in the read.
        Returns the tuple (insertions, deletions)
        insertions = [(pos1,ins1), (pos2,ins2)]
        posN = start position (preceding base, VCF-style)
        insN = length of inserted sequence (not including preceding base)
        deletions = [(pos1,del1), (pos2,del2)]
        posN = start position (preceding base, VCF-style)
        delN = length of deleted sequence (not including preceding base)
        Note: This does not count any "I" or "D" CIGAR operations at the start or end of a read.
        It also counts "N" as a deletion.
        """
        if self.__indels is None:
            blocks = self.get_contiguous_blocks( one_based=one_based )
            reverse = self.is_seq_reverse_complement()
            self.__indels = self._get_indels( blocks, reverse, one_based=one_based )
        return self.__indels

    @staticmethod
    def _get_indels( blocks, reverse, one_based=True ):
        #TODO: Include the cigar operation as a field in the block so this can avoid counting "N"s
        #      as deletions.
        insertions = []
        deletions = []
        last_read_end = None
        last_ref_end = None
        if reverse:
            for read_end, read_start, ref_end, ref_start, offset, direction in reversed(blocks):
                if last_read_end is not None:
                    if read_start == last_read_end:
                        del_len = last_ref_end-ref_start
                        del_start = last_ref_end-del_len-1
                        deletions.append( ( del_start, del_len ) )
                    else:
                        ins_start = last_ref_end-1
                        ins_len = read_start-last_read_end
                        insertions.append( ( ins_start, ins_len ) )
                last_read_end = read_end
                last_ref_end = ref_end
        else:
            for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
                if last_read_end is not None:
                    if read_start == last_read_end:
                        del_start = last_ref_end-1
                        del_len = ref_start-last_ref_end
                        deletions.append( ( del_start, del_len ) )
                    else:
                        ins_start = last_ref_end-1
                        ins_len = read_start-last_read_end
                        insertions.append( ( ins_start, ins_len ) )
                last_read_end = read_end
                last_ref_end = ref_end
        return (insertions, deletions)
    
    def get_end_position( self, one_based=True ):
        """Get the position of the end of the aligned portion of the sequence.
        This will always be the "right-most" position, regardless of orientation of the read.
        So for a read with its 5' end at 100 and its 3' end at 200, this will give 200.
        And for a reverse-oriented read with its 3' end at 100 and its 5' end at 200, this will
        still give 200.
        Note: hard- and soft-clipped bases are not counted."""
        if self.__zero_based_end_position is None:
            blocks = self.get_contiguous_blocks( one_based=False )
            self.__zero_based_end_position = self._get_end_position( blocks )
        return self.__zero_based_end_position + one_based

    @staticmethod
    def _get_end_position( blocks ):
        for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
            max_position = max(ref_start, ref_end)
        return max_position

    def get_5prime_position( self, one_based=True ):
        if self.is_seq_reverse_complement():
            return self.get_end_position( one_based=one_based )
        else:
            return self.get_position( one_based=True )

    def get_3prime_position( self, one_based=True ):
        if self.is_seq_reverse_complement():
            return self.get_position( one_based=True )
        else:
            return self.get_end_position( one_based=one_based )

    def get_pnext( self, one_based=True ):
        return self._next_pos + one_based
    def _get_bam_next_pos( self ):
        return pack_int32( self.get_pnext( False ) )
    
    def get_mapq( self ):
        return self._mapq
    
    def get_cigar( self ):
        self.__parse_block_2()
        return self._cigar_list[:]
    def get_sam_cigar( self ):
        self.__parse_block_2()
        return "".join( "%s%s" % ( l, CIGAR_OP[o] ) for l, o in self._cigar_list )
    
    def get_t_len( self ):
        return self._t_len
    def _get_bam_t_len( self ):
        return pack_int32( self._t_len )
    
    def get_seq( self ):
        if self.__seq_string is None:
            self.__parse_block_3()
            self.__seq_string = "".join( self._seq_list )
        return self.__seq_string
    def _get_bam_seq( self ):
        self.__parse_block_3()
        return struct.pack( "<" + "B" * (  ( self._l_seq + 1 ) / 2 ), *self._seq )

    def get_l_seq( self ):
        return self._l_seq
    def _get_bam_seq_length( self ):
        self.__parse_block_3()
        return pack_int32( len( self._seq_list ) )
    
    def get_qual( self ):
        self.__parse_block_4()
        return self._qual
    def get_qual_list( self ):
        if self.__qual_list is None:
            self.__qual_list = self.get_qual()
            if not self.__qual_list:
                self.__qual_list = [ BAM_NO_QUAL for i in range ( self._l_seq ) ]
        return self.__qual_list[:]
    def get_sam_qual( self ):
        return "".join( chr( c + 33 ) for c in self.get_qual() )
    def get_qual_tuple( self ):
        return tuple( self.get_qual_list() )
    def _get_bam_qual( self ):
        return struct.pack( "<" + "c" * self._l_seq, *map( chr, self.get_qual_tuple() ) )
    
    def _get_bam_bin_mq_nl( self ):
        return pack_uint32( self._get_bin_mq_nl() )
    
    def _get_bin_mq_nl( self ):
        return self._bin_mq_nl
    
    def _get_bam_flag_nc( self ):
        return pack_uint32( self._get_flag_nc() )
    def _get_flag_nc( self ):
        #TODO: fix me to calculate unless set and parsed already
        return self._flag_nc
    
    def _get_bam_n_cigar_op( self ):
        return struct.pack( "<" +"I" * self._n_cigar_op, *self._get_cigar() )
    
    def _get_cigar( self ):
        self.__parse_block_2()
        return self._cigar
    
    def is_seq_reverse_complement( self ):
        if self.__is_seq_reverse_complement is None:
            self.__is_seq_reverse_complement = ( self.get_flag() & 0x0010 == 0x0010 )
        return self.__is_seq_reverse_complement
    
    def get_read_group( self ):
        if self.__read_group_parsed is False:
            self.__parse_block_5()
            for tag, val_type, value in self._aux_data:
                if tag == READ_GROUP_RECORD_TAG:
                    self.__read_group = value[0]
                    break
            self.__read_group_parsed = True
        return self.__read_group
    def get_sam_aux( self ):
        self.__parse_block_5()
        rval = ''
        for aux in self._aux_data:
            rval += "\t%s:%s:%s" % ( aux[0], TAG_TYPE_TO_SAM_TYPE[ aux[1] ], ",".join( map( str, aux[2] ) ) )
        return rval.strip( '\t' )
    def _get_bam_aux( self ):
        self.__parse_block_5()
        data = ''
        for tag, val_type, value in self._aux_data:
            data += tag
            tag_length = len( value )
            if tag_length > 1:
                data = '%sB%s' ( data, val_type )
                data += pack_int32( tag_length )
            else:
                data += val_type
                if val_type in NULL_TERMINATED_TAGS:
                    data += value[0]
                    data += NULL_CHAR
                    tag_length = None
            if tag_length:
                data += struct.pack( "<" + TAG_TYPE_TO_STRUCT_TYPE[ val_type ] * tag_length, *value )
        return data
    
    def get_bam_data( self ):
        #FIX ME: have these calculate from updatable properties
        rval = "%s%s%s%s%s%s%s%s%s%s%s%s%s" % ( self._get_bam_ref_id(), self._get_bam_pos(), self._get_bam_bin_mq_nl(), 
                                                self._get_bam_flag_nc(), self._get_bam_seq_length(), self._get_bam_rnext_id(),
                                                self._get_bam_next_pos(), self._get_bam_t_len(), self._get_bam_read_name(),
                                                self._get_bam_n_cigar_op(), self._get_bam_seq(), self._get_bam_qual(),
                                                self._get_bam_aux() )
        rval = "%s%s" % ( pack_int32( len( rval ) ), rval )
        return rval
    
    def to_sam( self ):
        rval = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( self.get_read_name(), self.get_flag(), self.get_reference_name(), self.get_position(), self.get_mapq(), self.get_sam_cigar(), self.get_rnext_name(), self.get_pnext(), self.get_t_len(), self.get_seq(), self.get_sam_qual() )
        aux = self.get_sam_aux()
        if aux:
            rval = "%s\t%s" % ( rval, aux )
        return rval
