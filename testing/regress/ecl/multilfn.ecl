/*##############################################################################

    HPCC SYSTEMS software Copyright (C) 2012 HPCC Systems®.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
############################################################################## */

#onwarning(1005, ignore);

import std.system.thorlib;
import Std.File AS FileServices;
import std.str;

#option('slaveDaliClient', true);
#option('allFilenamesDynamic', true);

rec := RECORD, MAXLENGTH(256)
  STRING    album;
  STRING    song;
END;

ds := DATASET([
  {'BOB DYLAN', 'You\'re No Good'},
  {'BOB DYLAN', 'Talkin\' New York'},
  {'BOB DYLAN', 'In My Time Of Dyin\''},
  {'BOB DYLAN', 'Man Of Constant Sorrow'},
  {'BOB DYLAN', 'Fixin\' To Die'},
  {'BOB DYLAN', 'Pretty Peggy-O'},
  {'BOB DYLAN', 'Highway 51'},
  {'BOB DYLAN', 'Gospel Plow'},
  {'BOB DYLAN', 'Baby Let Me Follow You Down'},
  {'BOB DYLAN', 'House Of The Risin\' Sun'},
  {'BOB DYLAN', 'Freight Train Blues'},
  {'BOB DYLAN', 'Song To Woody'},
  {'BOB DYLAN', 'See That My Grave Is Kept Clean'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Blowin\' In The Wind'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Girl From The North Country'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Masters Of War'},
  {'THE FREEWHEELIN\' BOB DYLAN', ' Down The Highway'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Bob Dylan\'s Blues'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'A Hard Rain\'s A-Gonna Fall'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Don\'t Think Twice, It\'s All Right'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Bob Dylan\'s Dream'},
  {'THE FREEWHEELIN\' BOB DYLAN', ' Oxford Town'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Talking World War III Blues'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Corrina, Corrina'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'Honey, Just Allow Me One More Chance'},
  {'THE FREEWHEELIN\' BOB DYLAN', 'I Shall Be Free'},
  {'THE TIMES THEY ARE A-CHANGIN', 'The Times They Are A-Changin\''},
  {'THE TIMES THEY ARE A-CHANGIN', 'Ballad Of Hollis Brown'},
  {'THE TIMES THEY ARE A-CHANGIN', 'With God On Our Side'},
  {'THE TIMES THEY ARE A-CHANGIN', 'One Too Many Mornings'},
  {'THE TIMES THEY ARE A-CHANGIN', 'North Country Blues'},
  {'THE TIMES THEY ARE A-CHANGIN', 'Only A Pawn In Their Game'},
  {'THE TIMES THEY ARE A-CHANGIN', 'Boots Of Spanish Leather'},
  {'THE TIMES THEY ARE A-CHANGIN', 'When The Ship Comes In'},
  {'THE TIMES THEY ARE A-CHANGIN', 'Lonesome Death of Hattie Carroll'},
  {'THE TIMES THEY ARE A-CHANGIN', 'Restless Farewell'},
  {'ANOTHER SIDE OF BOB DYLAN', 'All I Really Want To Do'},
  {'ANOTHER SIDE OF BOB DYLAN', 'Black Crow Blues'},
  {'ANOTHER SIDE OF BOB DYLAN', 'Spanish Harlem Incident'},
  {'ANOTHER SIDE OF BOB DYLAN', 'Chimes Of Freedom'},
  {'ANOTHER SIDE OF BOB DYLAN', 'I Shall Be Free No.10'},
  {'ANOTHER SIDE OF BOB DYLAN', 'To Ramona'},
  {'ANOTHER SIDE OF BOB DYLAN', 'Motorpsycho Nitemare'},
  {'ANOTHER SIDE OF BOB DYLAN', 'My Back Pages'},
  {'ANOTHER SIDE OF BOB DYLAN', 'I Don\'t Believe You'},
  {'ANOTHER SIDE OF BOB DYLAN', 'Ballad In Plain D'},
  {'ANOTHER SIDE OF BOB DYLAN', 'It Ain\'t Me Babe'},
  {'BRINGING IT ALL BACK HOME', 'Subterranean Homesick Blues'},
  {'BRINGING IT ALL BACK HOME', 'She Belongs To Me'},
  {'BRINGING IT ALL BACK HOME', 'Maggie\'s Farm'},
  {'BRINGING IT ALL BACK HOME', 'Love Minus Zero/No Limit'},
  {'BRINGING IT ALL BACK HOME', 'Outlaw Blues'},
  {'BRINGING IT ALL BACK HOME', 'On The Road Again'},
  {'BRINGING IT ALL BACK HOME', 'Bob Dylan\'s 115th Dream'},
  {'BRINGING IT ALL BACK HOME', 'Mr. Tambourine Man'},
  {'BRINGING IT ALL BACK HOME', 'Gates Of Eden'},
  {'BRINGING IT ALL BACK HOME', 'It\'s Alright, Ma (I\'m Only Bleeding)'},
  {'BRINGING IT ALL BACK HOME', 'It\'s All Over Now, Baby Blue'},
  {'HIGHWAY 61 REVISITED', 'Like A Rolling Stone'},
  {'HIGHWAY 61 REVISITED', 'Tombstone Blues'},
  {'HIGHWAY 61 REVISITED', 'It Takes A Lot To Laugh, It Takes A Train To Cry'},
  {'HIGHWAY 61 REVISITED', 'From A Buick 6'},
  {'HIGHWAY 61 REVISITED', 'Ballad Of A Thin Man'},
  {'HIGHWAY 61 REVISITED', 'Queen Jane Approximately'},
  {'HIGHWAY 61 REVISITED', 'Highway 61 Revisited'},
  {'HIGHWAY 61 REVISITED', 'Just Like Tom Thumb\'s Blues'},
  {'HIGHWAY 61 REVISITED', 'Desolation Row'},
  {'BLONDE ON BLONDE', 'Rainy Day Women #12 & 35'},
  {'BLONDE ON BLONDE', 'Pledging My Time'},
  {'BLONDE ON BLONDE', 'Visions Of Johanna'},
  {'BLONDE ON BLONDE', 'One Of Us Must Know (Sooner Or Later)'},
  {'BLONDE ON BLONDE', 'I Want You'},
  {'BLONDE ON BLONDE', 'Stuck Inside Of Mobile With The Memphis Blues Again'},
  {'BLONDE ON BLONDE', 'Leopard-Skin Pill-Box Hat'},
  {'BLONDE ON BLONDE', 'Just Like A Woman'},
  {'BLONDE ON BLONDE', 'Most Likely You Go Your Way And I\'ll Go Mine'},
  {'BLONDE ON BLONDE', 'Temporary Like Achilles'},
  {'BLONDE ON BLONDE', 'Absolutely Sweet Marie'},
  {'BLONDE ON BLONDE', '4th Time Around'},
  {'BLONDE ON BLONDE', 'Obviously 5 Believers'},
  {'BLONDE ON BLONDE', 'Sad Eyed Lady Of The Lowlands'},
  {'JOHN WESLEY HARDING', 'John Wesley Harding'},
  {'JOHN WESLEY HARDING', 'As I Went Out One Morning'},
  {'JOHN WESLEY HARDING', 'I Dreamed I Saw St. Augustine'},
  {'JOHN WESLEY HARDING', 'All Along The Watchtower'},
  {'JOHN WESLEY HARDING', 'Ballad Of Frankie Lee & Judas Priest'},
  {'JOHN WESLEY HARDING', 'Drifter\'s Escape'},
  {'JOHN WESLEY HARDING', 'Dear Landlord'},
  {'JOHN WESLEY HARDING', 'I Am A Lonesome Hobo'},
  {'JOHN WESLEY HARDING', 'I Pity The Poor Immigrant'},
  {'JOHN WESLEY HARDING', 'The Wicked Messenger'},
  {'JOHN WESLEY HARDING', 'Down Along The Cove'},
  {'JOHN WESLEY HARDING', 'I\'ll Be Your Baby Tonight'},
  {'NASHVILLE SKYLINE', 'Girl From The North Country'},
  {'NASHVILLE SKYLINE', 'Nashville Skyline Rag'},
  {'NASHVILLE SKYLINE', 'To Be Alone With You'},
  {'NASHVILLE SKYLINE', 'I Threw It All Away'},
  {'NASHVILLE SKYLINE', 'Peggy Day'},
  {'NASHVILLE SKYLINE', 'Lay Lady Lay'},
  {'NASHVILLE SKYLINE', 'One More Night'},
  {'NASHVILLE SKYLINE', 'Tell Me That It Isn\'t True'},
  {'NASHVILLE SKYLINE', 'Country Pie'},
  {'NASHVILLE SKYLINE', 'Tonight I\'ll Be Staying Here With You'},
  {'SELF PORTRAIT', 'All The Tired Horses'},
  {'SELF PORTRAIT', 'Alberta #1'},
  {'SELF PORTRAIT', 'I Forgot More Than You\'ll Ever Know'},
  {'SELF PORTRAIT', 'Days Of 49'},
  {'SELF PORTRAIT', 'Early Morning Rain'},
  {'SELF PORTRAIT', 'In Search Of Little Sadie'},
  {'SELF PORTRAIT', 'Let It Be Me'},
  {'SELF PORTRAIT', 'Little Sadie'},
  {'SELF PORTRAIT', 'Woogie Boogie'},
  {'SELF PORTRAIT', 'Belle Isle'},
  {'SELF PORTRAIT', 'Living The Blues'},
  {'SELF PORTRAIT', 'Like A Rolling Stone'},
  {'SELF PORTRAIT', 'Copper Kettle'},
  {'SELF PORTRAIT', 'Gotta Travel On'},
  {'SELF PORTRAIT', 'Blue Moon'},
  {'SELF PORTRAIT', 'The Boxer'},
  {'SELF PORTRAIT', 'The Mighty Quinn'},
  {'SELF PORTRAIT', 'Take Me As I Am'},
  {'SELF PORTRAIT', 'Take A Message To Mary'},
  {'SELF PORTRAIT', 'It Hurts Me Too'},
  {'SELF PORTRAIT', 'Minstrel Boy'},
  {'SELF PORTRAIT', 'She Belongs To Me'},
  {'SELF PORTRAIT', 'Wigwam'},
  {'SELF PORTRAIT', 'Alberta #2'},
  {'NEW MORNING', 'If Not For You'},
  {'NEW MORNING', 'Day Of The Locusts'},
  {'NEW MORNING', 'Time Passes Slowly'},
  {'NEW MORNING', 'Went To See The Gypsy'},
  {'NEW MORNING', 'Winterlude'},
  {'NEW MORNING', 'If Dogs Run Free'},
  {'NEW MORNING', 'New Morning'},
  {'NEW MORNING', 'Sign On The Window'},
  {'NEW MORNING', 'One More Weekend'},
  {'NEW MORNING', 'The Man In Me'},
  {'NEW MORNING', 'Three Angels'},
  {'NEW MORNING', 'Father Of Night'},
  {'PAT GARRETT & BILLY THE KID', 'Main Title Theme (Billy)'},
  {'PAT GARRETT & BILLY THE KID', 'Cantina Theme (Workin\' For The Law)'},
  {'PAT GARRETT & BILLY THE KID', 'Billy 1'},
  {'PAT GARRETT & BILLY THE KID', 'Bunkhouse Theme'},
  {'PAT GARRETT & BILLY THE KID', 'River Theme'},
  {'PAT GARRETT & BILLY THE KID', 'Turkey Chase'},
  {'PAT GARRETT & BILLY THE KID', 'Knockin\' On Heaven\'s Door'},
  {'PAT GARRETT & BILLY THE KID', 'Final Theme'},
  {'PAT GARRETT & BILLY THE KID', 'Billy 4'},
  {'PAT GARRETT & BILLY THE KID', 'Billy 7'},
  {'DYLAN', 'Lily Of The West'},
  {'DYLAN', 'Can\'t Help Falling In Love'},
  {'DYLAN', 'Sarah Jane'},
  {'DYLAN', 'The Ballad Of Ira Hayes'},
  {'DYLAN', 'Mr. Bojangles'},
  {'DYLAN', 'Mary Ann'},
  {'DYLAN', 'Big Yellow Taxi'},
  {'DYLAN', 'A Fool Such As I'},
  {'DYLAN', 'Spanish Is The Loving Tongue'},
  {'PLANET WAVES', 'On A Night Like This'},
  {'PLANET WAVES', 'Going, Going, Gone'},
  {'PLANET WAVES', 'Tough Mama'},
  {'PLANET WAVES', 'Hazel'},
  {'PLANET WAVES', 'Something There Is About You'},
  {'PLANET WAVES', 'Forever Young (#1)'},
  {'PLANET WAVES', 'Forever Young (#2)'},
  {'PLANET WAVES', 'Dirge'},
  {'PLANET WAVES', 'You Angel You'},
  {'PLANET WAVES', 'Never Say Goodbye'},
  {'PLANET WAVES', 'Wedding Song'},
  {'BEFORE THE FLOOD', 'Most Likely You Go Your Way (And I\'ll Go Mine)'},
  {'BEFORE THE FLOOD', 'Lay Lady Lay'},
  {'BEFORE THE FLOOD', 'Rainy Day Women #12 & 35'},
  {'BEFORE THE FLOOD', 'Knockin\' On Heaven\'s Door'},
  {'BEFORE THE FLOOD', 'It Ain\'t Me Babe'},
  {'BEFORE THE FLOOD', 'Ballad Of A Thin Man'},
  {'BEFORE THE FLOOD', 'Up On Cripple Creek'},
  {'BEFORE THE FLOOD', 'I Shall Be Released'},
  {'BEFORE THE FLOOD', 'Endless Highway'},
  {'BEFORE THE FLOOD', 'The Night They Drove Old Dixie Down'},
  {'BEFORE THE FLOOD', 'Stage Fright'},
  {'BEFORE THE FLOOD', 'Don\'t Think Twice,It\'s All Right'},
  {'BEFORE THE FLOOD', 'Just Like A Woman'},
  {'BEFORE THE FLOOD', 'It\'s Alright, Ma (I\'m Only Bleeding)'},
  {'BEFORE THE FLOOD', 'The Shape I\'m In'},
  {'BEFORE THE FLOOD', 'When You Awake'},
  {'BEFORE THE FLOOD', 'The Weight'},
  {'BEFORE THE FLOOD', 'All Along The Watchtower'},
  {'BEFORE THE FLOOD', 'Highway 61 Revisited'},
  {'BEFORE THE FLOOD', 'Like A Rolling Stone'},
  {'BEFORE THE FLOOD', 'Blowin\' In The Wind'},
  {'BLOOD ON THE TRACKS', 'Tangled Up In Blue'},
  {'BLOOD ON THE TRACKS', 'Simple Twist Of Fate'},
  {'BLOOD ON THE TRACKS', 'You\'re A Big Girl Now'},
  {'BLOOD ON THE TRACKS', 'Idiot Wind'},
  {'BLOOD ON THE TRACKS', 'You\'re Gonna Make Me Lonesome When You Go'},
  {'BLOOD ON THE TRACKS', 'Meet Me In The Morning'},
  {'BLOOD ON THE TRACKS', 'Lily, Rosemary And The Of Hearts'},
  {'BLOOD ON THE TRACKS', 'If You See Her, Say Hello'},
  {'BLOOD ON THE TRACKS', 'Shelter From The Storm'},
  {'BLOOD ON THE TRACKS', 'Buckets Of Rain'},
  {'THE BASEMENT TAPES', 'Odds and Ends'},
  {'THE BASEMENT TAPES', 'Orange Juice Blues'},
  {'THE BASEMENT TAPES', 'Million Dollar Bash'},
  {'THE BASEMENT TAPES', 'Yazoo Street Scandal'},
  {'THE BASEMENT TAPES', 'Goin\' To Acapulco'},
  {'THE BASEMENT TAPES', 'Katie\'s Been Gone'},
  {'THE BASEMENT TAPES', 'Lo and Behold'},
  {'THE BASEMENT TAPES', 'Bessie Smith'},
  {'THE BASEMENT TAPES', 'Clothes Line Saga'},
  {'THE BASEMENT TAPES', 'Apple Suckling Tree'},
  {'THE BASEMENT TAPES', 'Please Mrs. Henry'},
  {'THE BASEMENT TAPES', 'Tears Of Rage'},
  {'THE BASEMENT TAPES', 'Too Much Of Nothing'},
  {'THE BASEMENT TAPES', 'Yea! Heavy and a Bottle of Bread'},
  {'THE BASEMENT TAPES', 'Ain\'t No More Cane'},
  {'THE BASEMENT TAPES', 'Crash On The Levee (Down In The Flood)'},
  {'THE BASEMENT TAPES', 'Ruben Remus'},
  {'THE BASEMENT TAPES', 'Tiny Montgomery'},
  {'THE BASEMENT TAPES', 'You Ain\'t Goin\' Nowhere'},
  {'THE BASEMENT TAPES', 'Don\'t Ya Tell Henry'},
  {'THE BASEMENT TAPES', 'Nothing Was Delivered'},
  {'THE BASEMENT TAPES', 'Open The Door, Homer'},
  {'THE BASEMENT TAPES', 'Long Distance Operator'},
  {'THE BASEMENT TAPES', 'This Wheel\'s On Fire'},
  {'DESIRE', 'Hurricane'},
  {'DESIRE', 'Isis'},
  {'DESIRE', 'Mozambique'},
  {'DESIRE', 'One More Cup Of Coffee'},
  {'DESIRE', 'Oh, Sister'},
  {'DESIRE', 'Joey'},
  {'DESIRE', 'Romance In Durango'},
  {'DESIRE', 'Black Diamond Bay'},
  {'DESIRE', 'Sara'},
  {'HARD RAIN', 'Maggie\'s Farm'},
  {'HARD RAIN', 'One Too Many Mornings'},
  {'HARD RAIN', 'Stuck Inside Of Mobile With The Memphis Blues Again'},
  {'HARD RAIN', 'Oh, Sister'},
  {'HARD RAIN', 'Lay Lady Lay'},
  {'HARD RAIN', 'Shelter From The Storm'},
  {'HARD RAIN', 'You\'re A Big Girl Now'},
  {'HARD RAIN', 'I Threw It All Away'},
  {'HARD RAIN', 'Idiot Wind'},
  {'STREET-LEGAL', 'Changing Of The Guards'},
  {'STREET-LEGAL', 'New Pony'},
  {'STREET-LEGAL', 'No Time To Think'},
  {'STREET-LEGAL', 'Baby Stop Crying'},
  {'STREET-LEGAL', 'Is Your Love In Vain?'},
  {'STREET-LEGAL', 'Se\361or (Tales Of Yankee Power)'},
  {'STREET-LEGAL', 'True Love Tends To Forget'},
  {'STREET-LEGAL', 'We Better Talk This Over'},
  {'STREET-LEGAL', 'Where Are You Tonight (Journey Through Dark Heat)'},
  {'AT BUDOKAN', 'Mr. Tambourine Man'},
  {'AT BUDOKAN', 'Shelter From The Storm'},
  {'AT BUDOKAN', 'Love Minus Zero/No Limit'},
  {'AT BUDOKAN', 'Ballad Of A Thin Man'},
  {'AT BUDOKAN', 'Don\'t Think Twice, It\'s All Right'},
  {'AT BUDOKAN', 'Maggie\'s Farm'},
  {'AT BUDOKAN', 'One More Cup Of Coffee (Valley Below)'},
  {'AT BUDOKAN', 'Like A Rolling Stone'},
  {'AT BUDOKAN', 'I Shall Be Released'},
  {'AT BUDOKAN', 'Is Your Love In Vain?'},
  {'AT BUDOKAN', 'Going, Going, Gone'},
  {'AT BUDOKAN', 'Blowin\' In The Wind'},
  {'AT BUDOKAN', 'Just Like A Woman'},
  {'AT BUDOKAN', 'Oh, Sister'},
  {'AT BUDOKAN', 'Simple Twist Of Fate'},
  {'AT BUDOKAN', 'All Along The Watchtower'},
  {'AT BUDOKAN', 'I Want You'},
  {'AT BUDOKAN', 'All I Really Want To Do'},
  {'AT BUDOKAN', 'Knockin\' On Heaven\'s Door'},
  {'AT BUDOKAN', 'It\'s Alright, Ma (I\'m Only Bleeding)'},
  {'AT BUDOKAN', 'Forever Young'},
  {'AT BUDOKAN', 'The Times They Are A-Changin\''},
  {'SLOW TRAIN COMING', 'Gotta Serve Somebody'},
  {'SLOW TRAIN COMING', 'Precious Angel'},
  {'SLOW TRAIN COMING', 'I Believe In You'},
  {'SLOW TRAIN COMING', 'Slow Train'},
  {'SLOW TRAIN COMING', 'Gonna Change My Way Of Thinking'},
  {'SLOW TRAIN COMING', 'Do Right To Me Baby'},
  {'SLOW TRAIN COMING', 'When You Gonna Wake Up'},
  {'SLOW TRAIN COMING', 'Man Gave Names To All The Animals'},
  {'SLOW TRAIN COMING', 'When He Returns'},
  {'SAVED', 'A Satisfied Mind'},
  {'SAVED', 'Saved'},
  {'SAVED', 'Covenant Woman'},
  {'SAVED', 'What Can I Do For You?'},
  {'SAVED', 'Solid Rock'},
  {'SAVED', 'Pressing On'},
  {'SAVED', 'In The Garden'},
  {'SAVED', 'Saving Grace'},
  {'SAVED', 'Are You Ready'},
  {'SHOT OF LOVE', 'Shot Of Love'},
  {'SHOT OF LOVE', 'Heart Of Mine'},
  {'SHOT OF LOVE', 'Property Of Jesus<'},
  {'SHOT OF LOVE', 'Lenny Bruce'},
  {'SHOT OF LOVE', 'Watered-Down Love'},
  {'SHOT OF LOVE', 'The Groom\'s Still Waiting At The Altar'},
  {'SHOT OF LOVE', 'Dead Man, Dead Man'},
  {'SHOT OF LOVE', 'In The Summertime'},
  {'SHOT OF LOVE', 'Trouble'},
  {'SHOT OF LOVE', 'Every Grain Of Sand'},
  {'INFIDELS', 'Jokerman'},
  {'INFIDELS', 'Sweetheart Like You'},
  {'INFIDELS', 'Neighborhood Bully'},
  {'INFIDELS', 'License To Kill'},
  {'INFIDELS', 'Man Of Peace'},
  {'INFIDELS', 'Union Sundown'},
  {'INFIDELS', 'I and I'},
  {'INFIDELS', 'Don\'t Fall Apart On Me Tonight'},
  {'REAL LIVE', 'Highway 61'},
  {'REAL LIVE', 'Maggie\'s Farm'},
  {'REAL LIVE', 'I and I'},
  {'REAL LIVE', 'License To Kill'},
  {'REAL LIVE', 'It Ain\'t Me Babe'},
  {'REAL LIVE', 'Tangled Up In Blue'},
  {'REAL LIVE', 'Masters Of War'},
  {'REAL LIVE', ' Ballad Of A Thin Man'},
  {'REAL LIVE', 'Girl From The North Country'},
  {'REAL LIVE', 'Tombstone Blues'},
  {'EMPIRE BURLESQUE', 'Tight Connection To My Heart'},
  {'EMPIRE BURLESQUE', 'Seeing The Real You At Last'},
  {'EMPIRE BURLESQUE', 'I\'ll Remember You'},
  {'EMPIRE BURLESQUE', 'Clean Cut Kid'},
  {'EMPIRE BURLESQUE', 'Never Gonna Be The Same Again'},
  {'EMPIRE BURLESQUE', 'Trust Yourself'},
  {'EMPIRE BURLESQUE', 'Emotionally Yours'},
  {'EMPIRE BURLESQUE', 'When The Night Comes Falling From The Sky'},
  {'EMPIRE BURLESQUE', 'Something\'s Burning, Baby'},
  {'EMPIRE BURLESQUE', 'Dark Eyes'},
  {'KNOCKED OUT LOADED', 'You Wanna Ramble'},
  {'KNOCKED OUT LOADED', 'They Killed Him'},
  {'KNOCKED OUT LOADED', 'Driftin\' Too Far From Shore'},
  {'KNOCKED OUT LOADED', 'Precious Memories'},
  {'KNOCKED OUT LOADED', 'Maybe Someday'},
  {'KNOCKED OUT LOADED', 'Brownsville Girl'},
  {'KNOCKED OUT LOADED', 'Got My Mind Made Up'},
  {'KNOCKED OUT LOADED', 'Under Your Spell'},
  {'DOWN IN THE GROOVE', 'Let\'s Stick Together'},
  {'DOWN IN THE GROOVE', 'When Did You Leave Heaven?'},
  {'DOWN IN THE GROOVE', 'Sally Sue Brown'},
  {'DOWN IN THE GROOVE', 'Death Is Not The End'},
  {'DOWN IN THE GROOVE', 'Had A Dream About You, Baby'},
  {'DOWN IN THE GROOVE', 'Ugliest Girl In The World'},
  {'DOWN IN THE GROOVE', 'Silvio'},
  {'DOWN IN THE GROOVE', 'Ninety Miles An Hour (Down A Dead End Street)'},
  {'DOWN IN THE GROOVE', 'Shenandoah'},
  {'DOWN IN THE GROOVE', 'Rank Strangers To Me'},
  {'DYLAN & THE DEAD', 'Slow Train'},
  {'DYLAN & THE DEAD', 'I Want You'},
  {'DYLAN & THE DEAD', 'Gotta Serve Somebody'},
  {'DYLAN & THE DEAD', 'Queen Jane Approximately'},
  {'DYLAN & THE DEAD', 'Joey'},
  {'DYLAN & THE DEAD', 'All Along The Watchtower'},
  {'DYLAN & THE DEAD', 'Knockin\' On Heaven\'s Door'},
  {'OH MERCY', 'Political World'},
  {'OH MERCY', 'Where Teardrops Fall'},
  {'OH MERCY', 'Everything Is Broken'},
  {'OH MERCY', 'Ring Them Bells'},
  {'OH MERCY', 'Man In The Long Black Coat'},
  {'OH MERCY', 'Most Of The Time'},
  {'OH MERCY', 'What Good Am I?'},
  {'OH MERCY', 'Disease Of Conceit'},
  {'OH MERCY', 'What Was It You Wanted'},
  {'OH MERCY', 'Shooting Star'},
  {'UNDER THE RED SKY', 'Wiggle Wiggle'},
  {'UNDER THE RED SKY', 'Under The Red Sky'},
  {'UNDER THE RED SKY', 'Unbelievable'},
  {'UNDER THE RED SKY', 'Born In Time'},
  {'UNDER THE RED SKY', 'T.V. Talkin\' Song'},
  {'UNDER THE RED SKY', '10,000 Men'},
  {'UNDER THE RED SKY', '2 X 2'},
  {'UNDER THE RED SKY', 'God Knows'},
  {'UNDER THE RED SKY', 'Handy Dandy'},
  {'UNDER THE RED SKY', 'Cat\'s In The Well'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Hard Times In New York Town'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'He Was A Friend Of Mine'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Man On The Street'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'No More Auction Block'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'House Carpenter'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Talkin\' Bear Mountain Picnic Massacre Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Let Me Die In My Footsteps'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Rambling, Gambling Willie'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Talkin\' Hava Negeilah Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Quit Your Low Down Ways'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Worried Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Kingsport Town'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Walkin\' Down The Line'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Walls of Red Wing'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Paths of Victory'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Talkin\' John Birch Paranoid Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Who Killed Davey Moore?'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Only A Hobo'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Moonshiner'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'When The Ship Comes In'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'The Times They Are A-Changin\''},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Last Thoughts on Woody Guthrie'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Seven Curses'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Eternal Circle'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Suze (The Cough Song)'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Mama, You Been On My Mind'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Farewell, Angelina'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Subterranean Homesick Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'If You Gotta Go, Go Now'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Sitting On A Barbed Wire Fence'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Like A Rolling Stone'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'It Takes A Lot To Laugh, It Takes A Train To Cry'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'I\'ll Keep It With Mine'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'She\'s Your Lover Now'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'I Shall Be Released'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Santa-Fe'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'If Not For You'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Wallflower'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Nobody \'Cept You'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Tangled Up In Blue'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Call Letter Blues'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Idiot Wind'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'If You See Her, Say Hello'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Golden Loom'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Catfish'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Seven Days'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Ye Shall Be Changed'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Every Grain Of Sand'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'You Changed My Life'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Need A Woman'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Angelina'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Someone\'s Got A Hold Of My Heart'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Tell Me'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Lord Protect My Child'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Foot Of Pride'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Blind Willie McTell'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'When The Night Comes Falling From The Sky'},
  {'THE BOOTLEG SERIES VOLUMES 1-3', 'Series Of Dreams'},
  {'GOOD AS I BEEN TO YOU', 'Frankie & Albert'},
  {'GOOD AS I BEEN TO YOU', 'Jim Jones'},
  {'GOOD AS I BEEN TO YOU', 'Blackjack Davey'},
  {'GOOD AS I BEEN TO YOU', 'Candee-I-O'},
  {'GOOD AS I BEEN TO YOU', 'Sittin\' On Top Of The World'},
  {'GOOD AS I BEEN TO YOU', 'Little Maggie'},
  {'GOOD AS I BEEN TO YOU', 'Hard Times'},
  {'GOOD AS I BEEN TO YOU', 'Step It Up And Go'},
  {'GOOD AS I BEEN TO YOU', 'Tomorrow Night'},
  {'GOOD AS I BEEN TO YOU', 'Arthur McBride'},
  {'GOOD AS I BEEN TO YOU', 'You\'re Gonna Quit Me'},
  {'GOOD AS I BEEN TO YOU', 'Diamond Joe'},
  {'GOOD AS I BEEN TO YOU', 'Froggie Went A Courtin\''},
  {'WORLD GONE WRONG', 'World Gone Wrong'},
  {'WORLD GONE WRONG', 'Love Henry'},
  {'WORLD GONE WRONG', 'Ragged & Dirty'},
  {'WORLD GONE WRONG', 'Blood In My Eyes'},
  {'WORLD GONE WRONG', 'Broke Down Engine'},
  {'WORLD GONE WRONG', 'Delia'},
  {'WORLD GONE WRONG', 'Stack A Lee'},
  {'WORLD GONE WRONG', 'Two Soldiers'},
  {'WORLD GONE WRONG', 'Jack-A-Roe'},
  {'WORLD GONE WRONG', 'Lone Pilgrim'},
  {'MTV UNPLUGGED', 'Tombstone Blues'},
  {'MTV UNPLUGGED', 'Shooting Star'},
  {'MTV UNPLUGGED', 'All Along The Watchtower'},
  {'MTV UNPLUGGED', 'The Times They Are A-Changin\''},
  {'MTV UNPLUGGED', 'John Brown'},
  {'MTV UNPLUGGED', 'Rainy Day Women #12 & 35'},
  {'MTV UNPLUGGED', 'Desolation Row'},
  {'MTV UNPLUGGED', 'Dignity'},
  {'MTV UNPLUGGED', 'Knockin\' On Heaven\'s Door'},
  {'MTV UNPLUGGED', 'Like A Rolling Stone'},
  {'MTV UNPLUGGED', 'With God On Our Side'},
  {'TIME OUT OF MIND', 'Love Sick'},
  {'TIME OUT OF MIND', 'Dirt Road Blues'},
  {'TIME OUT OF MIND', 'Standing In The Doorway'},
  {'TIME OUT OF MIND', 'Million Miles'},
  {'TIME OUT OF MIND', 'Tryin\' To Get To Heaven'},
  {'TIME OUT OF MIND', '\'Til I Fell In Love With You'},
  {'TIME OUT OF MIND', 'Not Dark Yet'},
  {'TIME OUT OF MIND', 'Cold Irons Bound'},
  {'TIME OUT OF MIND', 'Make You Feel My Love'},
  {'TIME OUT OF MIND', 'Can\'t Wait'},
  {'TIME OUT OF MIND', 'Highlands'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'She Belongs To Me'},
  {'THE "ROYAL ALBERT HALL" CONCERT', '4th Time Around'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Visions Of Johana'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'It\'s All Over Now, Baby Blue'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Desolation Row'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Mr. Tambourine Man'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Tell Me, Momma'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'I Don\'t Believe You'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Baby Let Me Follow You Down'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Just Like Tom Thumb\'s Blues'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Leopard-Skin Pill-Box Hat'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'One Too Many Mornings'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Ballad Of A Thin Man'},
  {'THE "ROYAL ALBERT HALL" CONCERT', 'Like A Rolling Stone'},
  {'LOVE AND THEFT', 'Tweedle Dee and Tweedle Dum'},
  {'LOVE AND THEFT', 'Mississippi'},
  {'LOVE AND THEFT', 'Summer Days'},
  {'LOVE AND THEFT', 'Bye and Bye'},
  {'LOVE AND THEFT', 'Lonesome Day Blues'},
  {'LOVE AND THEFT', 'Floater'},
  {'LOVE AND THEFT', 'Highwater (For Charley Patton)'},
  {'LOVE AND THEFT', 'Moonlight'},
  {'LOVE AND THEFT', 'Honest With Me'},
  {'LOVE AND THEFT', 'Po\' Boy'},
  {'LOVE AND THEFT', 'Cry Awhile'},
  {'LOVE AND THEFT', 'Sugar Baby'},
  {'THE ROLLING THUNDER REVUE', 'Tonight I\'ll Be Staying Here With You'},
  {'THE ROLLING THUNDER REVUE', 'It Ain\'t Me, Babe'},
  {'THE ROLLING THUNDER REVUE', 'A Hard Rain\'s A-Gonna Fall'},
  {'THE ROLLING THUNDER REVUE', 'The Lonesome Death Of Hattie Carroll'},
  {'THE ROLLING THUNDER REVUE', 'Romance In Durango'},
  {'THE ROLLING THUNDER REVUE', 'Isis'},
  {'THE ROLLING THUNDER REVUE', 'Mr. Tambourine Man'},
  {'THE ROLLING THUNDER REVUE', 'Simple Twist Of Fate'},
  {'THE ROLLING THUNDER REVUE', 'Blowin\' In The Wind'},
  {'THE ROLLING THUNDER REVUE', 'Mama, You Been On My Mind'},
  {'THE ROLLING THUNDER REVUE', 'I Shall Be Released'},
  {'THE ROLLING THUNDER REVUE', 'It\'s All Over Now, Baby Blue'},
  {'THE ROLLING THUNDER REVUE', 'Love Minus Zero/No Limit'},
  {'THE ROLLING THUNDER REVUE', 'Tangled Up In Blue'},
  {'THE ROLLING THUNDER REVUE', 'The Water Is Wide'},
  {'THE ROLLING THUNDER REVUE', 'It Takes A Lot To Laugh, It Takes A Train To Cry'},
  {'THE ROLLING THUNDER REVUE', 'Oh, Sister'},
  {'THE ROLLING THUNDER REVUE', 'Hurricane'},
  {'THE ROLLING THUNDER REVUE', 'One More Cup Of Coffee (Valley Below)'},
  {'THE ROLLING THUNDER REVUE', 'Sara'},
  {'THE ROLLING THUNDER REVUE', 'Just Like A Woman'},
  {'THE ROLLING THUNDER REVUE', 'Knockin\' On Heaven\'s Door'},
  {'CONCERT AT PHILHARMONIC HALL', 'The Times They Are A-Changin\''},
  {'CONCERT AT PHILHARMONIC HALL', 'Spanish Harlem Incident'},
  {'CONCERT AT PHILHARMONIC HALL', 'Talkin\' John Birch Paranoid Blues'},
  {'CONCERT AT PHILHARMONIC HALL', 'To Ramona'},
  {'CONCERT AT PHILHARMONIC HALL', 'Who Killed Davey Moore?'},
  {'CONCERT AT PHILHARMONIC HALL', 'Gates Of Eden'},
  {'CONCERT AT PHILHARMONIC HALL', 'If You Gotta Go, Go Now'},
  {'CONCERT AT PHILHARMONIC HALL', 'It\'s Alright, Ma (I\'m Only Bleeding)'},
  {'CONCERT AT PHILHARMONIC HALL', 'I Don\'t Believe You'},
  {'CONCERT AT PHILHARMONIC HALL', 'Mr. Tamborine Man'},
  {'CONCERT AT PHILHARMONIC HALL', 'A Hard Rain\'s A-Gonna Fall'},
  {'CONCERT AT PHILHARMONIC HALL', 'Talkin\' World War III Blues'},
  {'CONCERT AT PHILHARMONIC HALL', 'Don\'t Think Twice, It\'s All Right'},
  {'CONCERT AT PHILHARMONIC HALL', 'The Lonesome Death Of Hattie'},
  {'CONCERT AT PHILHARMONIC HALL', 'Mama, You Been On My Mind'},
  {'CONCERT AT PHILHARMONIC HALL', 'Silver Dagger'},
  {'CONCERT AT PHILHARMONIC HALL', 'With God On Our Side'},
  {'CONCERT AT PHILHARMONIC HALL', 'It Ain\'t Me, Babe'},
  {'CONCERT AT PHILHARMONIC HALL', 'All I Really Want To Do'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'When I Got Troubles'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Rambler, Gambler'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'This Land Is Your Land'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Song To Woody'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Dink\'s Song'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'I Was Young When I Left Home'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Sally Gal'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Don\'t Think Twice, It\'s All Right'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Man Of Constant Sorrow'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Blowin\' In The Wind'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Masters Of War'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'A Hard Rain\'s A-Gonna Fall'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'When The Ship Comes In'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Mr. Tambourine Man'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Chimes Of Freedom'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'It\'s All Over Now, Baby Blue'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'She Belongs To Me'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Maggie\'s Farm'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'It Takes A Lot To Laugh, It Takes A Train To Cry'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Tombstone Blues'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Just Like Tom Thumb\'s Blues'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Desolation Row'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Highway 61 Revisited'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Leopard-Skin Pill-Box Hat'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Stuck Inside Of Mobile With The Memphis Blues Again'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Visions Of Johanna'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Ballad Of A Thin Man'},
  {'NO DIRECTION HOME: THE SOUNDTRACK', 'Like A Rolling Stone'},
  {'MODERN TIMES', 'Thunder On The Mountain'},
  {'MODERN TIMES', 'Spirit on the Water'},
  {'MODERN TIMES', 'Rollin\' and Tumblin\''},
  {'MODERN TIMES', 'When the Deal Goes Down'},
  {'MODERN TIMES', 'Someday Baby'},
  {'MODERN TIMES', 'Workingman\'s Blues #2'},
  {'MODERN TIMES', 'Beyond the Horizon'},
  {'MODERN TIMES', 'Nettie Moore'},
  {'MODERN TIMES', 'The Levee\'s Gonna Break'},
  {'MODERN TIMES', 'Ain\'t Talkin\''}
], rec);


dds := DISTRIBUTE(ds,hash(album));
sds := SORT(dds,song);
o1 := output(sample(sds,5,1),,'A',overwrite);
o2 := output(sample(sds,5,2),,'~regress::BB',overwrite);
o3 := output(sample(sds,5,3),,'testscope::another::CCCC',overwrite);
o4 := output(sample(sds,5,4),,'~    regress :: DD DDD ',overwrite);   // space between DD and DDD significant
o5 := output(sample(sds,5,5),,'e^Ee^Ee',overwrite);

i1 := DATASET('{A,~regress::BB,testscope::another::cccc,~regress::dd ddd,e^Ee^Ee}',rec,flat);
o6 := OUTPUT(COUNT(ds)-COUNT(i1));
o7 := OUTPUT(COUNT(JOIN(i1,sds,left.song=right.song,FULL ONLY,LOCAL)));

o8 := output(sample(sds,5,1),,'nhtest_'+WORKUNIT+'::FFF',overwrite);
o9 := output(sample(sds,5,2),,'nhtest_'+WORKUNIT+'::GGG',overwrite);
o10 := output(sample(sds,5,3),,'nhtest_'+WORKUNIT+'::HHHHH',overwrite);
o11 := output(sample(sds,6,1),,'nhtest_'+WORKUNIT+'::HH',overwrite); // red herring
o12 := output(sample(sds,5,4),,'nhtest_'+WORKUNIT+'::III',overwrite);   // space between DD and DDD significant
o13 := output(sample(sds,5,5),,'nhtest_'+WORKUNIT+'::e^Ee^Ee',overwrite);

i2 := DATASET('nhtest_'+WORKUNIT+'::{???*}',rec,flat);
o14 := OUTPUT(COUNT(JOIN(i2,sds,left.song=right.song,FULL ONLY,LOCAL)));
o15 := OUTPUT(COUNT(i2)-COUNT(ds));

i3 := DATASET('~regress::{BB,dd ddd}',rec,flat);
o16 := OUTPUT(COUNT(JOIN(i3,sample(sds,5,2)+sample(sds,5,4),left.song=right.song,FULL ONLY)));

lfl1 := FileServices.LogicalFileList(str.ToLowerCase(thorlib.getExpandLogicalName('nhtest_'+WORKUNIT+'::*')));

a1 := NOTHOR(APPLY(lfl1,FileServices.DeleteLogicalFile('~'+name)));

SEQUENTIAL(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,a1);

