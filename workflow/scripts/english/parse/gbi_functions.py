import re

# eng book order
eng_book_list = '''
Genesis
Exodus
Leviticus
Numbers
Deuteronomy
Joshua
Judges
Ruth
1 Samuel
2 Samuel
1 Kings
2 Kings
1 Chronicles
2 Chronicles
Ezra
Nehemiah
Esther
Job
Psalms
Proverbs
Ecclesiastes
Song of songs
Isaiah
Jeremiah
Lamentations
Ezekiel
Daniel
Hosea
Joel
Amos
Obadiah
Jonah
Micah
Nahum
Habakkuk
Zephaniah
Haggai
Zechariah
Malachi
'''.strip().replace(' ', '_').split()

# map books to integers
int2book = {
    i+1: book for i, book in enumerate(eng_book_list)
}

# regex pattern for matching word ID info to its parts
# e.g. ('01', '001', '001', '001', '1')
# i.e. (bookN, chapterN, verseN, wordN, partN)
ref_id_re = re.compile('([0-9]{2})([0-9]{3})([0-9]{3})([0-9]{3})([1-9])')
ref_id_re_trans = re.compile('([0-9]{2})([0-9]{3})([0-9]{3})([0-9]{3})')

def id2ref(id_int, source='manuscript'):
    """Convert GBI ID ref tag to TF ref tuple"""
    
    id_str = str(id_int)
    
    # decide how to parse the ID tag
    methods = {
        'manuscript': {'re': ref_id_re, 'short': 11},
        'translation': {'re': ref_id_re_trans, 'short': 10}, 
    }
    
    method = methods[source]
    
    # fix ambiguity with lack of book padding in single-digit books
    if len(id_str) == method['short']:
        id_str = '0' + id_str
        
    # match reference parts
    try:
        bookN, chapterN, verseN = method['re'].match(id_str).groups()[:3]
    except AttributeError:
        raise Exception(f'wrong source for this id!')
    
    book = int2book.get(int(bookN))
    chapter = int(chapterN)
    verse = int(verseN)
    
    return (book, chapter, verse)
