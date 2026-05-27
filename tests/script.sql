-- Creer une table nono_body
--

CREATE TABLE nono_body (
  obs_id INTEGER NOT NULL,
  canal_id INTEGER NOT NULL,
  bt REAL ,
  o_a REAL,
  o_p REAL,
  mrq INTEGER ,
  E012233 REAL ,
  E055043 REAL ,
  E033032 INTEGER ,
  PRIMARY KEY(obs_id, canal_id),
  FOREIGN KEY(obs_id)
    REFERENCES airs_hdr(obs_id)
);

-- Creer une table nono_hdr
--

CREATE TABLE nono_hdr (
  obs_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  airs_id INTEGER NOT NULL,
  latitude REAL ,
  longitude REAL ,
  date INTEGER DEFAULT  20060101 ,
  time INTEGER DEFAULT  2300 ,
  orbite INTEGER ,
  delay INTEGER ,
  E005041 INTEGER, 
  E005043 INTEGER ,
  E007024 REAL ,
  E007025 REAL ,
  E055200 INTEGER ,
  E004209 INTEGER ,
  E005021 REAL DEFAULT 0.0,
  E005022 REAL ,
  E008012 INTEGER ,
  E025070 INTEGER ,
  E020010 INTEGER ,
  E007001 INTEGER ,
  FOREIGN KEY(airs_id)
    REFERENCES airs_parameter(airs_id)
);
