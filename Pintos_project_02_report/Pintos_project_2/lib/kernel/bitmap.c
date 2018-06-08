#include "bitmap.h"
#include <debug.h>
#include <limits.h>
#include <round.h>
#include <stdio.h>
#include "threads/malloc.h"
#ifdef FILESYS
#include "filesys/file.h"
#endif

/* Element type.

   This must be an unsigned integer type at least as wide as int.

   Each bit represents one bit in the bitmap.
   If bit 0 in an element represents bit K in the bitmap,
   then bit 1 in the element represents bit K+1 in the bitmap,
   and so on. */
typedef unsigned long elem_type;

/* Number of bits in an element. */
#define ELEM_BITS (sizeof (elem_type) * CHAR_BIT)

/* From the outside, a bitmap is an array of bits.  From the
   inside, it's an array of elem_type (defined above) that
   simulates an array of bits. */
/* (ADD) To contain bitmap index to indicate where to search next,
    add cur_index. This may reinitialize when allocating process
    is succressfully done, because it means some bits on bitmap
    was set, so index moves to bit next to allocated bits in last.*/
struct bitmap
  {
    size_t bit_cnt;     /* Number of bits. */
    elem_type *bits;    /* Elements that represent bits. */
	size_t cur_index;   /* (ADD) Index to start searching in Next fit strategy. */
  };

/* Returns the index of the element that contains the bit
   numbered BIT_IDX. */
static inline size_t
elem_idx (size_t bit_idx)
{
  return bit_idx / ELEM_BITS;
}

/* Returns an elem_type where only the bit corresponding to
   BIT_IDX is turned on. */
static inline elem_type
bit_mask (size_t bit_idx)
{
  return (elem_type) 1 << (bit_idx % ELEM_BITS);
}

/* Returns the number of elements required for BIT_CNT bits. */
static inline size_t
elem_cnt (size_t bit_cnt)
{
  return DIV_ROUND_UP (bit_cnt, ELEM_BITS);
}

/* Returns the number of bytes required for BIT_CNT bits. */
static inline size_t
byte_cnt (size_t bit_cnt)
{
  return sizeof (elem_type) * elem_cnt (bit_cnt);
}

/* Returns a bit mask in which the bits actually used in the last
   element of B's bits are set to 1 and the rest are set to 0. */
static inline elem_type
last_mask (const struct bitmap *b)
{
  int last_bits = b->bit_cnt % ELEM_BITS;
  return last_bits ? ((elem_type) 1 << last_bits) - 1 : (elem_type) -1;
}

/* Creation and destruction. */

/* Creates and returns a pointer to a newly allocated bitmap with room for
   BIT_CNT (or more) bits.  Returns a null pointer if memory allocation fails.
   The caller is responsible for freeing the bitmap, with bitmap_destroy(),
   when it is no longer needed. */
/* (ADD) Initialize cur_index = 0, means initial location to start searching. */
struct bitmap *
bitmap_create (size_t bit_cnt)
{
  struct bitmap *b = malloc (sizeof *b);
  if (b != NULL)
    {
      b->bit_cnt = bit_cnt;
	    b->cur_index = 0;    	  // Initialize current starting index to 0.
      b->bits = malloc (byte_cnt (bit_cnt));
      if (b->bits != NULL || bit_cnt == 0)
        {
          bitmap_set_all (b, false);
          return b;
        }
      free (b);
    }
  return NULL;
}

/* Creates and returns a bitmap with BIT_CNT bits in the
   BLOCK_SIZE bytes of storage preallocated at BLOCK.
   BLOCK_SIZE must be at least bitmap_needed_bytes(BIT_CNT). */
struct bitmap *
bitmap_create_in_buf (size_t bit_cnt, void *block, size_t block_size UNUSED)
{
  struct bitmap *b = block;

  ASSERT (block_size >= bitmap_buf_size (bit_cnt));

  b->bit_cnt = bit_cnt;
  b->bits = (elem_type *) (b + 1);
  bitmap_set_all (b, false);
  return b;
}

/* Returns the number of bytes required to accomodate a bitmap
   with BIT_CNT bits (for use with bitmap_create_in_buf()). */
size_t
bitmap_buf_size (size_t bit_cnt)
{
  return sizeof (struct bitmap) + byte_cnt (bit_cnt);
}

/* Destroys bitmap B, freeing its storage.
   Not for use on bitmaps created by bitmap_create_in_buf(). */
void
bitmap_destroy (struct bitmap *b)
{
  if (b != NULL)
    {
      free (b->bits);
      free (b);
    }
}

/* Bitmap size. */

/* Returns the number of bits in B. */
size_t
bitmap_size (const struct bitmap *b)
{
  return b->bit_cnt;
}

/* Setting and testing single bits. */

/* Atomically sets the bit numbered IDX in B to VALUE. */
void
bitmap_set (struct bitmap *b, size_t idx, bool value)
{
  ASSERT (b != NULL);
  ASSERT (idx < b->bit_cnt);
  if (value)
    bitmap_mark (b, idx);
  else
    bitmap_reset (b, idx);
}

/* Atomically sets the bit numbered BIT_IDX in B to true. */
void
bitmap_mark (struct bitmap *b, size_t bit_idx)
{
  size_t idx = elem_idx (bit_idx);
  elem_type mask = bit_mask (bit_idx);

  /* This is equivalent to `b->bits[idx] |= mask' except that it
     is guaranteed to be atomic on a uniprocessor machine.  See
     the description of the OR instruction in [IA32-v2b]. */
  asm ("orl %1, %0" : "=m" (b->bits[idx]) : "r" (mask) : "cc");
}

/* Atomically sets the bit numbered BIT_IDX in B to false. */
void
bitmap_reset (struct bitmap *b, size_t bit_idx)
{
  size_t idx = elem_idx (bit_idx);
  elem_type mask = bit_mask (bit_idx);

  /* This is equivalent to `b->bits[idx] &= ~mask' except that it
     is guaranteed to be atomic on a uniprocessor machine.  See
     the description of the AND instruction in [IA32-v2a]. */
  asm ("andl %1, %0" : "=m" (b->bits[idx]) : "r" (~mask) : "cc");
}

/* Atomically toggles the bit numbered IDX in B;
   that is, if it is true, makes it false,
   and if it is false, makes it true. */
void
bitmap_flip (struct bitmap *b, size_t bit_idx)
{
  size_t idx = elem_idx (bit_idx);
  elem_type mask = bit_mask (bit_idx);

  /* This is equivalent to `b->bits[idx] ^= mask' except that it
     is guaranteed to be atomic on a uniprocessor machine.  See
     the description of the XOR instruction in [IA32-v2b]. */
  asm ("xorl %1, %0" : "=m" (b->bits[idx]) : "r" (mask) : "cc");
}

/* Returns the value of the bit numbered IDX in B. */
bool
bitmap_test (const struct bitmap *b, size_t idx)
{
  ASSERT (b != NULL);
  ASSERT (idx < b->bit_cnt);
  return (b->bits[elem_idx (idx)] & bit_mask (idx)) != 0;
}

/* Setting and testing multiple bits. */

/* Sets all bits in B to VALUE. */
void
bitmap_set_all (struct bitmap *b, bool value)
{
  ASSERT (b != NULL);

  bitmap_set_multiple (b, 0, bitmap_size (b), value);
}

/* Sets the CNT bits starting at START in B to VALUE. */
void
bitmap_set_multiple (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t i;

  ASSERT (b != NULL);
  ASSERT (start <= b->bit_cnt);
  ASSERT (start + cnt <= b->bit_cnt);

  for (i = 0; i < cnt; i++)
    bitmap_set (b, start + i, value);
}

/* (ADD) To free frames allocated on buddy system, compute size that
   real amount of frame be set, and free them. */
void
bitmap_set_multiple_buddy (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t i;
  ASSERT (b != NULL);
  ASSERT (start <= b->bit_cnt);
  ASSERT (start + cnt <= b->bit_cnt);

  int j;
  int pref_size = 1;
  for (j = 256; j>= 1; j= j / 2)
  {
    if (cnt > j)
    {
      pref_size = j * 2;
      break;
    }
  }
  for (i = 0; i < pref_size; i++)
    bitmap_set (b, start + i, value);
}

/* Returns the number of bits in B between START and START + CNT,
   exclusive, that are set to VALUE. */
size_t
bitmap_count (const struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t i, value_cnt;

  ASSERT (b != NULL);
  ASSERT (start <= b->bit_cnt);
  ASSERT (start + cnt <= b->bit_cnt);

  value_cnt = 0;
  for (i = 0; i < cnt; i++)
    if (bitmap_test (b, start + i) == value)
      value_cnt++;
  return value_cnt;
}

/* Returns true if any bits in B between START and START + CNT,
   exclusive, are set to VALUE, and false otherwise. */
bool
bitmap_contains (const struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t i;

  ASSERT (b != NULL);
  ASSERT (start <= b->bit_cnt);
  ASSERT (start + cnt <= b->bit_cnt);

  for (i = 0; i < cnt; i++)
    if (bitmap_test (b, start + i) == value)
      return true;
  return false;
}

/* Returns true if any bits in B between START and START + CNT,
   exclusive, are set to true, and false otherwise.*/
bool
bitmap_any (const struct bitmap *b, size_t start, size_t cnt)
{
  return bitmap_contains (b, start, cnt, true);
}

/* Returns true if no bits in B between START and START + CNT,
   exclusive, are set to true, and false otherwise.*/
bool
bitmap_none (const struct bitmap *b, size_t start, size_t cnt)
{
  return !bitmap_contains (b, start, cnt, true);
}

/* Returns true if every bit in B between START and START + CNT,
   exclusive, is set to true, and false otherwise. */
bool
bitmap_all (const struct bitmap *b, size_t start, size_t cnt)
{
  return !bitmap_contains (b, start, cnt, false);
}

/* Finding set or unset bits. */

/* Finds and returns the starting index of the first group of CNT
   consecutive bits in B at or after START that are all set to
   VALUE.
   If there is no such group, returns BITMAP_ERROR. */
/* (ADD) Scan fuction for First fit. */
/* ALLOCATOR_FF, 0.*/
size_t
bitmap_scan (const struct bitmap *b, size_t start, size_t cnt, bool value)
{
  ASSERT (b != NULL);
  ASSERT (start <= b->bit_cnt);

  if (cnt <= b->bit_cnt)
    {
      size_t last = b->bit_cnt - cnt;
      size_t i;
      for (i = start; i <= last; i++)
        if (!bitmap_contains (b, i, cnt, !value))
          return i;
    }
  return BITMAP_ERROR;
}
/* ALLOCATOR_NF, 1.*/
/* (ADD) Finds and returns the starting index of next group of CNT
	bit. To search from start index if searching index reached to last,
  with no finding bits to be allocated, iterator moves starting index
  and search again when it reached start index initially. */
size_t
bitmap_scan_nf(struct bitmap *b, size_t start, size_t cnt, bool value)
{
	ASSERT(b != NULL);
	ASSERT(start <= b->bit_cnt);

	if (cnt <= b->bit_cnt)
	{
		size_t last = b->bit_cnt - cnt;
		size_t i;
		start = b->cur_index;
		for (i = start; i <= last; i++)
			if (!bitmap_contains(b, i, cnt, !value))
			{
				b->cur_index = i+cnt;
				return i;
			}
		for(i = 0; i < start - cnt; i++){
			if (!bitmap_contains(b, i, cnt, !value)) 
			{
				b->cur_index = i + cnt;
				return i;
			}
      }
	}
	return BITMAP_ERROR;
}

/* ALLOCATOR_BF, 2.*/
/* (ADD) Find and return the starting index that has most profit size.
   To search smallest profit size, it saves index and increases current size,
   until it meets set bit. When it finds set bit, compare current size with
   min_size (current smallest size ever), if current size is smaller than
   min_size, min_size become current size and current start index becomes
   min_index. To do this process, initial min_size declared on 512, never be
   smaller than any size, min_index on 0 to allocate when all bits in bitmap
   is 0. */
size_t
bitmap_scan_bf(const struct bitmap *b, size_t start, size_t cnt, bool value)
{
	ASSERT(b != NULL);
	ASSERT(start <= b->bit_cnt);
	
	if (cnt <= b->bit_cnt)
	{
		size_t last = b->bit_cnt - cnt;
		size_t i;
		size_t j;
		size_t min_size = 512;
		size_t min_index = start;
		size_t size = 0;
		for (i = start; i <= last; i++)
		{
			if (!bitmap_contains(b, i, cnt, !value))
			{
				for (j = i; j < b->bit_cnt && !bitmap_test(b, j); j++)
				{
					size++;
				}
				if (size < min_size)
				{
					min_index = i;
					min_size = size;
					size = 0;
				}
				i = j;
			}
		}
		return min_index;
	}
	return BITMAP_ERROR;
}

/* ALLOCATOR_BUDDY, 3.*/
/* (ADD) Finds starting index to allocate and return that. Especially, in
   in buddy system, allocating unit is limited power of 2, and its size is
   computed 2^i <= size needed < 2^(i+1). Upon 'I', found in that fomular,
   buddy system allocates 2^(i+1) size, iteratively dividing so get prefer size,
   and that process would be shown as tree shape with 2 leap node. But in this
   allocating algorithm, no such tree struct is needed. Cause this system
   allocate more frames that process really needs, free in same way,
   it would be sufficient to search bitmap by iterator that has prefer size.
   To free bigger frames that allocated, buddy system needs other free procedure,
   implemented in method bitmap_set_multiple_buddy(), and it must be called in
   page_free function when allocate mode is buddy system. */
size_t
bitmap_scan_buddy(const struct bitmap *b, size_t start, size_t cnt, bool value)
{
	ASSERT(b != NULL);
	ASSERT(start <= b->bit_cnt);

	if (cnt <= b->bit_cnt)
	{
		int j;
		int pref_size = 1;
		for (j = 256; j>= 1; j= j / 2)
		{
			if (cnt > j)
			{
				pref_size = j * 2;
        break;
			}
		}

		size_t last = b->bit_cnt - pref_size;
		size_t i;
		for (i = start; i <= last; i += pref_size)
			if (!bitmap_contains(b, i, pref_size, !value))
				return i;
	}
	return BITMAP_ERROR;
}

/* Finds the first group of CNT consecutive bits in B at or after
   START that are all set to VALUE, flips them all to !VALUE,
   and returns the index of the first bit in the group.
   If there is no such group, returns BITMAP_ERROR.
   If CNT is zero, returns 0.
   Bits are set atomically, but testing bits is not atomic with
   setting them. */
size_t
bitmap_scan_and_flip (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t idx = bitmap_scan (b, start, cnt, value);
  if (idx != BITMAP_ERROR)
    bitmap_set_multiple (b, idx, cnt, !value);
  return idx;
}

/* (ADD) Same fucntion in First fit, just scan on Next fit way. */
size_t
bitmap_scan_and_flip_nf (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t idx = bitmap_scan_nf (b, start, cnt, value);
  if (idx != BITMAP_ERROR)
    bitmap_set_multiple (b, idx, cnt, !value);
  return idx;
}

/* (ADD) Same fucntion in First fit, just scan on Best fit way. */
size_t
bitmap_scan_and_flip_bf (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  size_t idx = bitmap_scan_bf (b, start, cnt, value);
  if (idx != BITMAP_ERROR)
    bitmap_set_multiple (b, idx, cnt, !value);
  return idx;
}
/* (ADD) To set all bits may allocated for current process, compute prefer size
    and set all of them. Size of free frames in guaranteed to be fit because
    bitmap_scan_buddy() also compute same value. */
size_t
bitmap_scan_and_flip_buddy (struct bitmap *b, size_t start, size_t cnt, bool value)
{
  int j;
  int pref_size = 1;
  size_t idx = bitmap_scan_buddy (b, start, cnt, value);
  if (idx != BITMAP_ERROR)
  {
    for (j = 256; j>= 1; j= j / 2)
		{
			if (cnt > j)
			{
				pref_size = j * 2;
        break;
			}
		}
    bitmap_set_multiple (b, idx, pref_size, !value);
  }
  return idx;
}

/* File input and output. */

#ifdef FILESYS
/* Returns the number of bytes needed to store B in a file. */
size_t
bitmap_file_size (const struct bitmap *b)
{
  return byte_cnt (b->bit_cnt);
}

/* Reads B from FILE.  Returns true if successful, false
   otherwise. */
bool
bitmap_read (struct bitmap *b, struct file *file)
{
  bool success = true;
  if (b->bit_cnt > 0)
    {
      off_t size = byte_cnt (b->bit_cnt);
      success = file_read_at (file, b->bits, size, 0) == size;
      b->bits[elem_cnt (b->bit_cnt) - 1] &= last_mask (b);
    }
  return success;
}

/* Writes B to FILE.  Return true if successful, false
   otherwise. */
bool
bitmap_write (const struct bitmap *b, struct file *file)
{
  off_t size = byte_cnt (b->bit_cnt);
  return file_write_at (file, b->bits, size, 0) == size;
}
#endif /* FILESYS */

/* Debugging. */

/* Dumps the contents of B to the console as hexadecimal. */
void
bitmap_dump (const struct bitmap *b)
{
  hex_dump (0, b->bits, byte_cnt (b->bit_cnt), false);
}

/* Dumps the contents of B to the console as binary. */
void
bitmap_dump2 (const struct bitmap *b)
{
  size_t i, j;

  for (i=0; i<elem_cnt (b->bit_cnt); i++) {
    for (j=0; j<ELEM_BITS; j++) {
      if ((i * ELEM_BITS + j) < b->bit_cnt) {
        printf ("%u", (unsigned int) (b->bits[i] >> j) & 0x1);
      }
    }
    printf ("\n");
  }
}
